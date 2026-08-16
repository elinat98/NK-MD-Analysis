# MD Trajectory Analysis — Hydrogen Bonds, Representative Frames, RMSD

Analysis notebooks for multi-replicate molecular dynamics trajectories of a
protein–protein interface: hydrogen bond occupancy per residue, representative
frame selection, and Cα RMSD.

| Notebook | Purpose | Input | Key outputs |
|---|---|---|---|
| `01_hbond_analysis.ipynb` | Per-residue and total interface H-bonds between two chain groups, both chemical directions | PDB + DCDs | 3 CSVs, 4 figures, `.pkl` results object |
| `02_representative_frame_analysis.ipynb` | Frame best representing the ensemble | `.pkl` from 01 | Representative frame PDB, ranked alternatives |
| `03_rmsd_analysis.ipynb` | Per-frame Cα RMSD, time series and distributions | PDB + DCDs | RMSD CSVs and figures |

Notebook 02 depends on notebook 01. Notebook 03 is standalone and uses the same
path configuration format.

---

## Requirements

Python 3.10+, MDAnalysis, NumPy, pandas, Matplotlib, seaborn. MDAnalysis and
seaborn are installed in-notebook via `pip`.

Written for **Google Colab** with trajectories on Google Drive. To run locally,
delete the `drive.mount(...)` call and the `files.download(...)` calls and point
the paths at a local directory — nothing else changes.

---

## Configuration

Replace the `INSERT_...` placeholders at the top of each notebook:

```python
pdb_file = 'INSERT_TOPOLOGY_FILENAME.pdb'
dcd_files = [
    'INSERT_TRAJECTORY_FILENAME_R1.dcd',
    'INSERT_TRAJECTORY_FILENAME_R2.dcd',
    'INSERT_TRAJECTORY_FILENAME_R3.dcd',
]
output_dir = 'INSERT_OUTPUT_DIRECTORY/'
```

The interface is defined as two chain groups:

```python
CHAINS_1 = ['C']
TARGET_RESIDUES_1 = []        # [] -> whole chain group, no per-residue breakdown

CHAINS_2 = ['A', 'B']
TARGET_RESIDUES_2 = [147, 148, ...]
```

The search is symmetric — each group donates to and accepts from the other — so
the groups are interchangeable. Group 1 is listed first in residue-pair labels
and appears on the left of the frequency plots; that is the only difference.

Each group's `TARGET_RESIDUES_` list is optional:

- **Non-empty** — the *target* scope is restricted to those residues, and each
  gets its own bar in the per-residue figures.
- **Empty (`[]`)** — the whole chain group takes part in the target scope, and
  no per-residue breakdown is produced for that side.

With both lists empty the target scope is identical to the total scope: plain
H-bonding between the two chain groups. Residue numbers are resids as they
appear in the PDB topology.

---

## Methods as implemented

### Hydrogen bond criteria

Maximum donor–acceptor distance **3.5 Å**, minimum donor–hydrogen–acceptor angle
**150°**.

Donor, hydrogen and acceptor atoms are given as explicit atom-name lists rather
than left to automatic guessing, so definitions are identical across replicates
and directions.

### Search direction

Both chemical directions are searched — group 1 → group 2 and group 2 → group 1
— since each group contains donors and acceptors. Two scopes are analysed per
replicate, giving four searches in total:

- **target** — the residues selected on each side by the `TARGET_RESIDUES_` lists
- **total** — the two chain groups in full, whole interface

### Counting

Reported frequencies are **atom-level**: a residue pair forming two hydrogen
bonds in one frame contributes 2. Units are therefore *bonds per frame* and
values can exceed 1. A single bond increments the count of the participating
residue on each side that has a target list.

Separately, each frame's residue-level pairs are stored as a set — a pair counts
once per frame however many atom-level bonds it forms, in either direction.
These sets are used **only** by notebook 02 for representative-frame scoring, and
are keyed `(group 1 residue, group 2 residue)` whichever side donated, so the
two searches describe the same pair rather than mirror-image entries.

### Equilibration

The first **100 ns** of each replicate is discarded. The frame count is derived
from a time so the value stays correct if the save interval changes:

```python
FRAME_INTERVAL_NS = 0.1       # coordinates saved every 100 ps
EQUILIBRATION_NS = 100.0      # 100 ns -> FRAMES_TO_SKIP = 1000
```

With coordinates saved every 100 ps, a 1 µs replicate is 10,000 frames and
100 ns is 1,000 of them. Across three replicates that leaves **27,000 of 30,000
frames** analysed.

The skip is written into the results object, so notebook 02 inherits it
(`FRAMES_TO_SKIP = None`) and its occupancies match the frequency tables from
notebook 01 by construction. Set the same `EQUILIBRATION_NS` in notebook 03.

**Where the skip is and is not applied.** Per-frame arrays — `frame_pairs`,
`rmsd_values`, `hbond_counts` — are recorded for the *full* trajectory with
absolute frame indices; the skip is applied when statistics are computed, not
when data is recorded. A frame number therefore means the same thing in
notebook 01, notebook 02, and any PDB extracted from a trajectory.

`RMSD_AVERAGE_START` separately controls the RMSD reference: `0` (default)
averages over the full trajectory including equilibration; set it to
`FRAMES_TO_SKIP` to average over the production portion only.

### Representative frame selection

For each frame *i* after the equilibration skip:

```
characteristic = { pair : freq  |  freq >= MIN_FREQ }        # MIN_FREQ = 5%

identity_i  = sum(freq of characteristic pairs present in frame i)
              / sum(freq of all characteristic pairs)

rmsd_norm_i = min-max normalized Ca RMSD over the pooled scored frames

score_i     = IDENTITY_WEIGHT * (1 - identity_i) + RMSD_WEIGHT * rmsd_norm_i
```

The lowest-scoring frame wins — maximum frequency-weighted identity, minimum
structural deviation. Defaults: `MIN_FREQ = 5.0`, `IDENTITY_WEIGHT = 0.7`,
`RMSD_WEIGHT = 0.3`. Scoring is direction-agnostic. The top four frames are
printed so an alternative can be exported by hand.

### RMSD

`calculate_rmsd_to_average` computes the Cα RMSD of each frame to its own
replicate's average structure, using a Kabsch/SVD superposition. Notebook 03
also offers `RMSD_MODE = 'reference'`, which measures against one fixed frame
using the identical superposition, so the two traces are on the same scale.
