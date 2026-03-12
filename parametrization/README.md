# Parametrization Resources

This folder collects the parametrization material used for the AGE crosslink examples in this repository.

## Contents

### [`Antechamber_steps`](./Antechamber_steps/)

This folder contains the step-by-step Antechamber workflow and the worked glucosepane example.

- [`Antechamber_steps/README.md`](./Antechamber_steps/README.md): commands and notes for the Antechamber parametrization.
- [`Antechamber_steps/glucosepane/`](./Antechamber_steps/glucosepane/): generated files for the glucosepane example.

### [`topology_antechamber+Grappa`](./topology_antechamber+Grappa/)

This folder contains ready-to-use topology files for each crosslink chemistry using two different parametrization schemes:

- `topol.top`: topology obtained with the Antechamber-based parametrization.
- `grappa_topol.top`: topology obtained with Grappa starting from the Antechamber topology.

The folder currently includes paired topologies for:

- `glucosepane`
- `pentosidine`
- `MOLD`

It also includes the updated `amber99sb-star-ildnp.ff.zip` force-field directory with the additional parameters needed to run MD with the Antechamber-derived topologies in the Amber99sb*-ildnp framework.

- [`topology_antechamber+Grappa/README.md`](./topology_antechamber+Grappa/README.md): details on the topology files and how they were generated.
