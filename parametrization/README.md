# Glucosepane Parametrization

This folder documents how we parametrized the glucosepane crosslink.

**Note:** The workflow assumes you already ran a QM calculation for glucosepane. We used Gaussian16 to generate the required output file before starting the steps below. The Gaussian input file with the required keywords and associated IOPs is provided as`glucosepane.com`.

## Dependencies

Install the following locally before starting:

1. **AmberTools (latest)**

```bash
conda create --name AmberTools25 python=3.12
conda activate AmberTools25
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install dacase::ambertools-dac=25
```

2. **acpype**

```bash
conda install -c conda-forge acpype
```

3. **Gaussian16**

Run a QM calculation for glucosepane and make sure the output file is available (see the first command below).

## Parametrization Steps

```bash
espgen -i glucosepane.out -o glucosepane.esp
```

```bash
antechamber -i glucosepane.out -fi gout -o glucosepane.ac -fo ac
```

```bash
respgen -i glucosepane.ac -o glucosepane.respin1 -f resp1
```

```bash
respgen -i glucosepane.ac -o glucosepane.respin2 -f resp2
```

```bash
emacs glucosepane.respin1
emacs glucosepane.respin2
emacs glucosepane.qin
```

Edit these three files as follows:
- In `glucosepane.respin1` and `glucosepane.respin2`, set `ivary = -1` for the atoms belonging to the capping groups you want to freeze.
- In `glucosepane.qin`, enter the frozen charges for those capping-group atoms at the corresponding positions.
  The capping-group charges are listed in `capping_charges`.

```bash
resp -O -i glucosepane.respin1 -o glucosepane.respout1 -e glucosepane.esp -t qout_stage1 -q glucosepane.qin
```

```bash
resp -O -i glucosepane.respin2 -o glucosepane.respout2 -e glucosepane.esp -q qout_stage1 -t qout_stage2
```

```bash
antechamber -i glucosepane.ac -fi ac -o glucosepane_resp.ac -fo ac -c rc -cf qout_stage2
```

```bash
atomtype -i glucosepane_resp.ac -o glucosepane_resp_gaff.ac -f ac -p gaff
```

```bash
prepgen -i glucosepane_resp_gaff.ac -o glucosepane.prepi -f car
```

```bash
parmchk -i glucosepane.prepi -f prepc -o glucosepane.frcmod
```

```bash
antechamber -i glucosepane.prepi -fi prepc -o glucosepane_resp_gaff.mol2 -fo mol2
```

```bash
acpype -p glucosepane.prmtop -x glucosepane.inpcrd -g
```

After `acpype`, the `.gro`, `.top` and `.itp` files will be written to `MOL.amb2gmx` as:
- `MOL_GMX.gro`
- `MOL_GMX.top`
- `posre_MOL.itp`
