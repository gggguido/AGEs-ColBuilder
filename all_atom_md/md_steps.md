# MD Workflow Steps

This guide describes how to run MD simulations on collagen microfibrils produced by ColBuilder. The topology generation step in ColBuilder produces the following files, all of which are used in the workflow below:

- `col_*.itp`
- `posre_*.itp`
- `collagen_fibril_rattus_norvegicus.gro`
- `collagen_fibril_rattus_norvegicus.top`
- Force field directory: `amber99sb-star-ildnp.ff`

## System Preparation

1. Center and align fibril

```bash
gmx editconf -f collagen_fibril_rattus_norvegicus.gro -o fibril_x.gro -c -princ <<<1
```

2. Rotate

```bash
gmx editconf -f fibril_x.gro -o fibril_z.gro -rotate 0 270 0 -c
```

3. Define box

```bash
gmx editconf -f fibril_z.gro -o fibril_box.gro -c -box 19 20 95 -bt triclinic
```

4. Solvate

```bash
gmx solvate -cp fibril_box.gro -p collagen_fibril_rattus_norvegicus.top -o fibril_solv.gro
```

5. Generate ions input
- Required mdp file location: `minimization/ions.mdp`

```bash
gmx grompp -f minimization/ions.mdp -c fibril_solv.gro -p collagen_fibril_rattus_norvegicus.top -o fibril_genion.tpr
```

6. Add ions
- Output written to: `minimization/ions.gro`

```bash
gmx genion -s fibril_genion.tpr -p collagen_fibril_rattus_norvegicus.top -o minimization/ions.gro -conc 0.15 -neutral < SOL
```

7. Vacuum minimization (setup + run)
- Required mdp file location: `minimization/minim_vacuum.mdp`

```bash
gmx grompp -f minimization/minim_vacuum.mdp -c minimization/ions.gro -p collagen_fibril_rattus_norvegicus.top -o fibril_min_vacuum.tpr
gmx mdrun -deffnm fibril_min_vacuum -v -maxwarn 1
```

8. Solvated minimization (setup + run)
- Required mdp file location: `minimization/minim_solv.mdp`

```bash
gmx grompp -f minimization/minim_solv.mdp -c fibril_min_vacuum.gro -p collagen_fibril_rattus_norvegicus.top -o fibril_min_solv.tpr
gmx mdrun -deffnm fibril_min_solv -v
```

## Equilibration

**Position restraints note (posre_*.itp):**

ColBuilder generates `posre_*.itp` files with a default force constant of **1.0 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>**. In this workflow, we used different position-restraint force constants at specific stages, so the `posre_*.itp` files were updated accordingly:

- **A.NVT, 0.5 fs timestep:** 2.5 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>
  - Use `update_posre_2500.py` in `Equilibration/A.NVT/2500/` to update the force constant.
- **B.NPT, 1 fs and 2 fs timesteps:** 0.3 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>
  - Use `update_posre_300.py` in `Equilibration/B.NPT/300/` to update the force constant.

Make sure to adjust the scripts' `EXPECTED_FILES` values to match the number of `posre_*.itp` files in your system.


### NVT (stepwise timestep and restraints)

mdp locations:
- `Equilibration/A.NVT/2500/0.5dt/nvt.mdp`
- `Equilibration/A.NVT/1000/1dt/2nvt.mdp`
- `Equilibration/A.NVT/1000/2dt/3nvt.mdp`

**A) 0.5 fs timestep, position restraint force constant 2.5 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>**

```bash
gmx grompp -f nvt.mdp -c fibril_min_solv.gro -r fibril_min_solv.gro -p collagen_fibril_rattus_norvegicus.top -o fibril_nvt.tpr
gmx mdrun -deffnm fibril_nvt -v -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1
```

**B) 1.0 fs timestep, position restraint force constant 1.0 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>**

```bash
gmx grompp -f 2nvt.mdp -c fibril_nvt.gro -r fibril_nvt.gro -p collagen_fibril_rattus_norvegicus.top -t fibril_nvt.cpt -o fibril_nvt.tpr
gmx mdrun -deffnm fibril_nvt -v -cpi fibril_nvt.cpt -append -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1
```

**C) 2.0 fs timestep, position restraint force constant 1.0 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>**

```bash
gmx grompp -f 3nvt.mdp -c fibril_nvt.gro -r fibril_nvt.gro -p collagen_fibril_rattus_norvegicus.top -t fibril_nvt.cpt -o fibril_nvt.tpr
gmx mdrun -deffnm fibril_nvt -v -cpi fibril_nvt.cpt -append -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1
```

### NPT (stepwise timestep and restraints)

mdp locations:
- `Equilibration/B.NPT/1000/0.5dt/npt.mdp`
- `Equilibration/B.NPT/300/1dt/2npt.mdp`
- `Equilibration/B.NPT/300/2dt/3npt.mdp`

**A) 0.5 fs timestep, position restraint force constant 1.0 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>**

```bash
gmx grompp -f npt.mdp -c fibril_nvt.gro -r fibril_nvt.gro -t fibril_nvt.cpt -p collagen_fibril_rattus_norvegicus.top -o fibril_npt.tpr
gmx mdrun -deffnm fibril_npt -v -noappend -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1
```

**B) 1.0 fs timestep, position restraint force constant 0.3 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>**

```bash
gmx grompp -f 2npt.mdp -c fibril_npt.part0001.gro -r fibril_npt.part0001.gro -t fibril_npt.cpt -p collagen_fibril_rattus_norvegicus.top -o fibril_npt.tpr
gmx mdrun -deffnm fibril_npt.part0001 -s fibril_npt.tpr -cpi fibril_npt.cpt -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1
```

**C) 2.0 fs timestep, position restraint force constant 0.3 kcal·mol<sup>−1</sup>·Å<sup>−2</sup>**

```bash
gmx grompp -f 3npt.mdp -c fibril_npt.part0002.gro -r fibril_npt.part0002.gro -t fibril_npt.part0001.cpt -p collagen_fibril_rattus_norvegicus.top -o fibril_npt.tpr
gmx mdrun -deffnm fibril_npt.part0002 -s fibril_npt.tpr -cpi fibril_npt.part0001.cpt -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1
```


## Production Run

mdp files, example `index.ndx`, and the in-house script (`enforced_rotation.py`) for setting up enforced rotation groups are available in `Production`.

### Index and Groups

1. Create index file

```bash
gmx make_ndx -f fibril_npt.part0003.gro -o index.ndx
```

2. Define groups

- Select `1|13|14` and name the group `all_dry`.
  - `1` = Protein, `13` = AGS, `14` = LGX.
  - For glucosepane, AGS and LGX are the markers that define the divalent glucosepane.
  - General logic: include all crosslink markers present in the collagen microfibril.
- Select `r ACE & a CH3` and name the group `ACE_&_CH3`.
- Select `r NME & a CH3` and name the group `NME_&_CH3`.

3. Generate appended index and mdp fragments

```bash
python enforced_rotation.py index.ndx
```

This script reads `index.ndx` and expects the `ACE_&_CH3` and `NME_&_CH3` groups. It splits each group into chunks of three atoms, creates per-chunk rotation groups named `ACE_0`, `ACE_1`, ... and `NME_0`, `NME_1`, ..., then writes:

- `append.ndx` with the new rotation groups (to be appended to `index.ndx`).
- `append.mdp` with enforced-rotation settings (`rotation = Yes`, `rot-ngroups = <count>`, `rot-type = rm-pf`, `rot-vec = 0 0 1`, `rot-rate = 0.0`, `rot-k = 200`, `rot-fit-method = norm`) for those groups.

Append `append.ndx` to `index.ndx` and append `append.mdp` to the production mdp file.

### Pulling Protocol (5 steps, +0.2 nN each to 1.0 nN)

We used a gradual force ramp to reach 1.0 nN in five steps. This is optional; you can apply the target force in step 1 and run only the first step if you prefer.

**Step 1**

```bash
gmx grompp -f md-pull_step1.mdp -c fibril_npt.part0003.gro -t fibril_npt.part0002.cpt -o fibril_pull.tpr -p collagen_fibril_rattus_norvegicus.top -n index.ndx -maxwarn 2
gmx mdrun -deffnm fibril_pull -v -noappend -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1
```

**Step 2**

```bash
gmx grompp -f md-pull_step2.mdp -c fibril_pull.part0001.gro -t fibril_pull.cpt -o fibril_pull_part2.tpr -p collagen_fibril_rattus_norvegicus.top -n index.ndx -maxwarn 1
gmx mdrun -s fibril_pull_part2.tpr -deffnm fibril_pull -cpi fibril_pull.cpt -cpo fibril_pull.cpt -noappend -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1 -v
```

**Step 3**

```bash
gmx grompp -f md-pull_step3.mdp -c fibril_pull.part0002.gro -t fibril_pull.cpt -o fibril_pull_part3.tpr -p collagen_fibril_rattus_norvegicus.top -n index.ndx -maxwarn 1
gmx mdrun -s fibril_pull_part3.tpr -deffnm fibril_pull -cpi fibril_pull.cpt -cpo fibril_pull.cpt -noappend -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1 -v
```

**Step 4**

```bash
gmx grompp -f md-pull_step4.mdp -c fibril_pull.part0003.gro -t fibril_pull.cpt -o fibril_pull_part4.tpr -p collagen_fibril_rattus_norvegicus.top -n index.ndx -maxwarn 1
gmx mdrun -s fibril_pull_step4.tpr -deffnm fibril_pull -cpi fibril_pull.cpt -cpo fibril_pull.cpt -noappend -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1 -v
```

**Step 5**

```bash
gmx grompp -f md-pull_step5.mdp -c fibril_pull.part0004.gro -t fibril_pull.cpt -o fibril_pull_part5.tpr -p collagen_fibril_rattus_norvegicus.top -n index.ndx -maxwarn 1
gmx mdrun -s fibril_pull_part5.tpr -deffnm fibril_pull -cpi fibril_pull.cpt -cpo fibril_pull.cpt -noappend -ntmpi 4 -ntomp 5 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -npme 1 -v
```

Before each step, repeat steps 2 and 3 using the `.gro` file produced at the end of the previous step to recreate `index.ndx` and regenerate the rotation groups.
