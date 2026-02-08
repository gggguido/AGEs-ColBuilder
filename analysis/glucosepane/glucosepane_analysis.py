#!/usr/bin/env python3
"""
Analyze collagen microfibril extension, overlap, and gap from ACE/NME end caps.

Method summary:
- Per frame, compute end-to-end distance as the Euclidean distance between the
  mean positions of ACE (N-terminus) and NME (C-terminus); this is the extension
  (D-band length proxy).
- Identify crosslink marker sets from the bond network (AGS/LGX pairs (glucosepane)) and compute
  each marker-set mean position per frame.
- Project marker-set positions onto the instantaneous ACE→NME axis and apply
  k=2 1D K-means to split N- and C-terminal crosslink clusters.
- Define overlap as the distance between cluster centroids; define gap as
  end-to-end minus overlap.
- Estimate overlap/gap uncertainties from within-cluster spread (SEM) and propagate
  to overlap/gap; for end-to-end, combine propagated uncertainty with a sliding
  window bootstrap (default window=10, resamples=1000), combined in quadrature.
- Compute overlap/gap strain ratio relative to the first analyzed frame and
  propagate uncertainty.
- Average across replicas and combine per-replica errors in quadrature.

Helper (duo markers via m1/m2):
  python glucosepane_analysis.py \
  --top reduced.tpr \
  --traj traj_vmd.xtc \
  --sel1 "resname ACE" --sel2 "resname NME" \
  --m1-resname AGS --m2-resname LGX \
  --start 0 --stop -1 --stride 1 \
  --bootstrap-window 10 --bootstrap-resamples 1000 \
  --out-prefix collagen67nm
"""
import argparse
import csv
import sys
from collections import defaultdict, deque

import numpy as np
import MDAnalysis as mda

AGS_RESNAME = "AGS"
LGX_RESNAME = "LGX"


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Clustering-based (k=2) analysis across replicas: mean positions, SEM from within-cluster spread, error propagation, and sliding-window bootstrap. "
            "Crosslink = AGS+LGX duo identified via bonds."
        )
    )
    p.add_argument("--top", required=True, help="Topology consistent with XTC (TPR recommended: contains bonds)")
    p.add_argument("--traj", required=True, nargs="+", help="XTC trajectories (e.g., 3 independent replicas)")

    p.add_argument("--sel1", default="resname ACE", help="N-terminal group (default: resname ACE)")
    p.add_argument("--sel2", default="resname NME", help="C-terminal group (default: resname NME)")

    p.add_argument("--m1-resname", default=AGS_RESNAME, help="Marker 1 resname (default: AGS)")
    p.add_argument("--m2-resname", default=LGX_RESNAME, help="Marker 2 resname (default: LGX)")

    p.add_argument("--start", type=int, default=0, help="First frame (inclusive)")
    p.add_argument("--stop", type=int, default=-1, help="Last frame (exclusive), -1 = end")
    p.add_argument("--stride", type=int, default=1, help="Analyze every N frames (1 = all frames)")
    p.add_argument("--dt-ps", type=float, default=None, help="Fallback dt [ps] if ts.time is missing")

    # Bootstrap parameters
    p.add_argument("--bootstrap-window", type=int, default=10, help="Bootstrap window (frames), default=10")
    p.add_argument("--bootstrap-resamples", type=int, default=1000, help="Resamples per window, default=1000")
    p.add_argument("--bootstrap-seed", type=int, default=0, help="Bootstrap RNG seed")

    p.add_argument("--out-prefix", default="clustering_k2", help="Output prefix")

    return p.parse_args()


def get_time_ps(ts, u, dt_ps_fallback):
    if getattr(ts, "time", None) is not None:
        return float(ts.time)
    dt = getattr(u.trajectory, "dt", None)
    if dt is not None:
        return float(ts.frame) * float(dt)
    if dt_ps_fallback is not None:
        return float(ts.frame) * float(dt_ps_fallback)
    return float(ts.frame)


def connected_components(nodes, adj):
    seen = set()
    comps = []
    for n in nodes:
        if n in seen:
            continue
        q = deque([n])
        seen.add(n)
        comp = {n}
        while q:
            x = q.popleft()
            for y in adj.get(x, ()):
                if y in nodes and y not in seen:
                    seen.add(y)
                    comp.add(y)
                    q.append(y)
        comps.append(comp)
    return comps


def build_marker_residue_graph(u, marker_resnames):
    if not hasattr(u, "bonds") or u.bonds is None or len(u.bonds) == 0:
        raise ValueError("Topology has no bonds (u.bonds empty). A bonded topology is required (e.g., TPR).")

    marker_set = set(marker_resnames)
    marker_nodes = set(r.ix for r in u.residues if r.resname in marker_set)

    adj = defaultdict(set)
    for b in u.bonds:
        a1, a2 = b.atoms[0], b.atoms[1]
        r1, r2 = a1.residue, a2.residue
        if r1 is r2:
            continue
        if (r1.ix in marker_nodes) and (r2.ix in marker_nodes):
            adj[r1.ix].add(r2.ix)
            adj[r2.ix].add(r1.ix)

    return adj, marker_nodes


def detect_duos(u, m1_resname="AGS", m2_resname="LGX"):
    """
    A valid "duo" is a connected component (using ONLY bonds among marker residues)
    containing exactly 2 residues: 1x m1 and 1x m2.
    """
    adj, marker_nodes = build_marker_residue_graph(u, [m1_resname, m2_resname])
    comps = connected_components(marker_nodes, adj)

    duos = []
    bad = 0
    for comp in comps:
        residues = [u.residues[ix] for ix in comp]
        m1_list = [r for r in residues if r.resname == m1_resname]
        m2_list = [r for r in residues if r.resname == m2_resname]

        if len(residues) == 2 and len(m1_list) == 1 and len(m2_list) == 1:
            duos.append((m1_list[0], m2_list[0]))
        else:
            bad += 1

    duos.sort(key=lambda d: (int(d[0].resid), int(d[1].resid)))
    return duos, bad, len(comps)


def kmeans_1d_k2(x, n_iter=50):
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        return np.zeros_like(x, dtype=int)

    c0 = float(np.min(x))
    c1 = float(np.max(x))
    if c0 == c1:
        return np.zeros_like(x, dtype=int)

    for _ in range(n_iter):
        d0 = np.abs(x - c0)
        d1 = np.abs(x - c1)
        lab = (d1 < d0).astype(int)

        if np.all(lab == 0) or np.all(lab == 1):
            c0, c1 = float(np.min(x)), float(np.max(x))
            if c0 == c1:
                return np.zeros_like(x, dtype=int)
            continue

        new_c0 = float(x[lab == 0].mean())
        new_c1 = float(x[lab == 1].mean())

        if (abs(new_c0 - c0) < 1e-12) and (abs(new_c1 - c1) < 1e-12):
            break

        c0, c1 = new_c0, new_c1

    return lab


def sem_of_mean_1d(values):
    v = np.asarray(values, dtype=float)
    n = int(v.size)
    if n <= 1:
        return np.nan
    return float(np.std(v, ddof=1) / np.sqrt(n))


def sem_of_center_xyz(coords_nm):
    n = int(coords_nm.shape[0])
    if n <= 1:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    std_xyz = np.std(coords_nm, axis=0, ddof=1)
    return std_xyz / np.sqrt(n)


def sem_distance_from_centers(p1, sem1_xyz, p2, sem2_xyz):
    r = p2 - p1
    d = float(np.linalg.norm(r))
    if d == 0.0:
        return np.nan
    var_r = sem1_xyz**2 + sem2_xyz**2
    grad = r / d
    var_d = float(np.sum((grad**2) * var_r))
    return float(np.sqrt(var_d))


def bootstrap_sem_sliding(values, window, B, seed):
    rng = np.random.default_rng(seed)
    x = np.asarray(values, dtype=float)
    T = int(x.size)
    out = np.full(T, np.nan, dtype=float)

    if window < 2 or B < 2 or T < window:
        return out

    for i in range(window - 1, T):
        w = x[i - window + 1 : i + 1]
        idx = rng.integers(0, window, size=(B, window))
        means = w[idx].mean(axis=1)
        out[i] = float(np.std(means, ddof=1))

    return out


def run_one_traj(top, traj, sel1, sel2, m1_resname, m2_resname, start, stop, stride, dt_ps_fallback):
    u = mda.Universe(top, traj)

    cap1 = u.select_atoms(sel1)
    cap2 = u.select_atoms(sel2)
    if len(cap1) == 0 or len(cap2) == 0:
        raise ValueError(f"Empty selections: len(sel1)={len(cap1)} len(sel2)={len(cap2)}")

    duos, bad, ncomps = detect_duos(u, m1_resname, m2_resname)
    if len(duos) == 0:
        raise ValueError("No valid duos found (check marker resnames and bonded topology).")
    if bad > 0:
        print(f"[WARN] {traj}: discarded anomalous marker components: {bad} (tot={ncomps})", file=sys.stderr)

    n_duos = len(duos)
    if n_duos % 2 != 0:
        raise ValueError(f"Odd number of duos: {n_duos} (must be even for a balanced k=2 split).")

    duo_idx = []
    for (r1, r2) in duos:
        idx = np.unique(np.concatenate([r1.atoms.indices, r2.atoms.indices]))
        duo_idx.append(idx)

    n_frames = len(u.trajectory)
    i0 = max(0, start)
    i1 = n_frames if stop is None or stop < 0 else min(stop, n_frames)
    st = max(1, stride)
    if i0 >= i1:
        raise ValueError(f"Invalid frame range: start={i0} stop={i1}")

    times = []
    all_nm = []
    overlap_nm = []
    gap_nm = []

    overlap_sem = []
    all_sem_geo = []
    gap_sem = []

    u.trajectory[i0]
    for ts in u.trajectory[i0:i1:st]:
        t_ns = get_time_ps(ts, u, dt_ps_fallback) * 1e-3
        coords_nm = u.atoms.positions * 0.1  # Angstrom -> nm

        c1 = coords_nm[cap1.indices]
        c2 = coords_nm[cap2.indices]

        p1 = c1.mean(axis=0)  # mean position ACE
        p2 = c2.mean(axis=0)  # mean position NME

        r = p2 - p1
        d_all = float(np.linalg.norm(r))
        if d_all == 0.0:
            continue
        axis = r / d_all

        sem1_xyz = sem_of_center_xyz(c1)
        sem2_xyz = sem_of_center_xyz(c2)
        s_all_geo = sem_distance_from_centers(p1, sem1_xyz, p2, sem2_xyz)

        # 1D projections of duos along ACE->NME axis
        s = np.empty(n_duos, dtype=float)
        for i in range(n_duos):
            pd = coords_nm[duo_idx[i]].mean(axis=0)  # mean position duo
            s[i] = float(np.dot(pd - p1, axis))

        lab = kmeans_1d_k2(s, n_iter=50)

        mu0 = float(np.mean(s[lab == 0]))
        mu1 = float(np.mean(s[lab == 1]))
        if mu0 <= mu1:
            sN = s[lab == 0]
            sC = s[lab == 1]
        else:
            sN = s[lab == 1]
            sC = s[lab == 0]

        muN = float(np.mean(sN))
        muC = float(np.mean(sC))
        d_ov = float(abs(muC - muN))

        semN = sem_of_mean_1d(sN)
        semC = sem_of_mean_1d(sC)
        s_ov = float(np.sqrt(semN**2 + semC**2))

        d_gap = float(d_all - d_ov)
        s_gap = float(np.sqrt((s_all_geo**2) + (s_ov**2)))

        times.append(float(t_ns))
        all_nm.append(d_all)
        overlap_nm.append(d_ov)
        gap_nm.append(d_gap)

        all_sem_geo.append(s_all_geo)
        overlap_sem.append(s_ov)
        gap_sem.append(s_gap)

    return {
        "times": np.asarray(times, float),
        "all": np.asarray(all_nm, float),
        "overlap": np.asarray(overlap_nm, float),
        "gap": np.asarray(gap_nm, float),
        "all_sem_geo": np.asarray(all_sem_geo, float),
        "overlap_sem": np.asarray(overlap_sem, float),
        "gap_sem": np.asarray(gap_sem, float),
    }


def write_mean_file(path, t, mean, err):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_ns", "mean", "error"])
        for ti, m, e in zip(t, mean, err):
            w.writerow([f"{ti:.6f}", f"{m:.6f}", f"{e:.6f}"])


def main():
    args = parse_args()
    n_runs = len(args.traj)

    runs = []
    for tr in args.traj:
        print(f"[INFO] Analyzing: {tr}", file=sys.stderr)
        r = run_one_traj(
            top=args.top,
            traj=tr,
            sel1=args.sel1,
            sel2=args.sel2,
            m1_resname=args.m1_resname,
            m2_resname=args.m2_resname,
            start=args.start,
            stop=args.stop,
            stride=args.stride,
            dt_ps_fallback=args.dt_ps,
        )
        runs.append(r)

    L = min(len(r["times"]) for r in runs)
    if L == 0:
        print("[ERROR] No data produced.", file=sys.stderr)
        sys.exit(2)

    t = runs[0]["times"][:L]

    A = np.vstack([r["all"][:L] for r in runs])
    O = np.vstack([r["overlap"][:L] for r in runs])
    G = np.vstack([r["gap"][:L] for r in runs])

    A_geo = np.vstack([r["all_sem_geo"][:L] for r in runs])
    O_sem = np.vstack([r["overlap_sem"][:L] for r in runs])
    G_sem = np.vstack([r["gap_sem"][:L] for r in runs])

    # D-band / all error per run: (A) propagated overlap+gap + (B) bootstrap, RMS
    A_prop = np.sqrt(O_sem**2 + G_sem**2)
    A_boot = []
    for i in range(n_runs):
        b = bootstrap_sem_sliding(
            values=A[i, :],
            window=args.bootstrap_window,
            B=args.bootstrap_resamples,
            seed=args.bootstrap_seed + i,
        )
        A_boot.append(b)
    A_boot = np.vstack(A_boot)
    A_tot = np.sqrt(A_prop**2 + A_boot**2)

    # Per-run ratio: (O(t)-O0)/(G(t)-G0) and propagated uncertainty
    R = np.full_like(O, np.nan, dtype=float)
    R_sem = np.full_like(O, np.nan, dtype=float)

    for i in range(n_runs):
        O0 = float(O[i, 0])
        G0 = float(G[i, 0])
        sO0 = float(O_sem[i, 0])
        sG0 = float(G_sem[i, 0])

        dO = O[i, :] - O0
        dG = G[i, :] - G0

        s_dO = np.sqrt(O_sem[i, :]**2 + sO0**2)
        s_dG = np.sqrt(G_sem[i, :]**2 + sG0**2)

        with np.errstate(divide="ignore", invalid="ignore"):
            R[i, :] = dO / dG
            R_sem[i, :] = np.abs(R[i, :]) * np.sqrt((s_dO / dO) ** 2 + (s_dG / dG) ** 2)

        R[i, ~np.isfinite(R[i, :])] = np.nan
        R_sem[i, ~np.isfinite(R_sem[i, :])] = np.nan

    # Means across runs
    A_mean = np.nanmean(A, axis=0)
    O_mean = np.nanmean(O, axis=0)
    G_mean = np.nanmean(G, axis=0)
    R_mean = np.nanmean(R, axis=0)

    # Error on the replica mean (independent replicas):
    # sigma_mean = sqrt(sum_i sigma_i^2) / n_runs
    def err_on_mean(S):
        return np.sqrt(np.nansum(S**2, axis=0)) / float(n_runs)

    A_err = err_on_mean(A_tot)
    O_err = err_on_mean(O_sem)
    G_err = err_on_mean(G_sem)
    R_err = err_on_mean(R_sem)

    prefix = args.out_prefix
    write_mean_file(f"{prefix}_mean_all.csv", t, A_mean, A_err)
    write_mean_file(f"{prefix}_mean_overlap.csv", t, O_mean, O_err)
    write_mean_file(f"{prefix}_mean_gap.csv", t, G_mean, G_err)
    write_mean_file(f"{prefix}_mean_ratio.csv", t, R_mean, R_err)

    print("[OK] wrote:", file=sys.stderr)
    print(f"  {prefix}_mean_all.csv", file=sys.stderr)
    print(f"  {prefix}_mean_overlap.csv", file=sys.stderr)
    print(f"  {prefix}_mean_gap.csv", file=sys.stderr)
    print(f"  {prefix}_mean_ratio.csv", file=sys.stderr)

    # Debug per-run
    with open(f"{prefix}_debug_run_errors.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "run", "time_ns",
            "all_nm", "all_sem_geo", "all_sem_prop", "all_sem_boot", "all_sem_total",
            "overlap_nm", "overlap_sem",
            "gap_nm", "gap_sem",
            "ratio", "ratio_sem"
        ])
        for irun in range(n_runs):
            for k in range(L):
                w.writerow([
                    irun + 1,
                    f"{t[k]:.6f}",
                    f"{A[irun, k]:.6f}",
                    f"{A_geo[irun, k]:.6f}",
                    f"{A_prop[irun, k]:.6f}",
                    f"{A_boot[irun, k]:.6f}",
                    f"{A_tot[irun, k]:.6f}",
                    f"{O[irun, k]:.6f}",
                    f"{O_sem[irun, k]:.6f}",
                    f"{G[irun, k]:.6f}",
                    f"{G_sem[irun, k]:.6f}",
                    f"{R[irun, k]:.6f}" if np.isfinite(R[irun, k]) else "",
                    f"{R_sem[irun, k]:.6f}" if np.isfinite(R_sem[irun, k]) else "",
                ])


if __name__ == "__main__":
    main()
