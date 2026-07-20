# Structure-Based Analysis of Flavonoid Selectivity for Chronic Wound Healing Targets Using Molecular Docking and Interaction Analysis

*In simpler words: analyzing flavonoids to find which ones could potentially help in chronic wound healing.*

---

## Table of Contents

- [Abstract](#abstract)
- [Background and Rationale](#background-and-rationale)
- [Objectives](#objectives)
- [Protein Targets](#protein-targets)
- [Tools and Software Used](#tools-and-software-used)
- [Workflow](#workflow)
  - [1. Environment Setup](#1-environment-setup)
  - [2. Ligand Acquisition](#2-ligand-acquisition)
  - [3. Ligand Merging](#3-ligand-merging)
  - [4. Ligand Renaming](#4-ligand-renaming)
  - [5. Ligand Preparation](#5-ligand-preparation)
  - [6. Protein Structure Acquisition](#6-protein-structure-acquisition)
  - [7. Protein Cleaning](#7-protein-cleaning)
  - [8. Reference Ligand Isolation and Processing](#8-reference-ligand-isolation-and-processing)
  - [9. Protein Preparation for Docking](#9-protein-preparation-for-docking)
  - [10. Grid Box Validation](#10-grid-box-validation)
  - [11. Docking Setup](#11-docking-setup)
  - [12. Molecular Docking](#12-molecular-docking)
  - [13. Extracting Binding Affinities](#13-extracting-binding-affinities)
  - [14. Converting Poses to SDF](#14-converting-poses-to-sdf)
  - [15. Pose Selection for Interaction Analysis](#15-pose-selection-for-interaction-analysis)
  - [16. Protein-Ligand Interaction Analysis](#16-protein-ligand-interaction-analysis)
  - [17. Selectivity Analysis](#17-selectivity-analysis)
- [Docking Results](#docking-results)
- [Interaction Analysis](#interaction-analysis)
  - [COX-2 (5IKR)](#cox-2-5ikr)
  - [MMP-9 (6ESM)](#mmp-9-6esm)
  - [EGFR (6VHN)](#egfr-6vhn)
  - [VEGFR-2 (3VO3)](#vegfr-2-3vo3)
- [Comparative Discussion](#comparative-discussion)
  - [Isoquercetin](#isoquercetin)
  - [Kaempferol](#kaempferol)
- [Limitations](#limitations)
- [Proposed Next Steps](#proposed-next-steps)

---

## Abstract

This project aims to discover potential drug candidates among flavonoids for chronic wound healing by using molecular docking and interaction analysis to find compounds that could potentially **selectively inhibit** inflammatory and ECM-degrading receptors while **sparing** vital wound-healing receptors involved in epithelial cell growth and angiogenesis. This is done by docking flavonoids and their glycones against **MMP-9**, **COX-2**, **EGFR**, and **VEGFR-2**. The docked poses are then analyzed for the residues they interact with, in order to predict their inhibitory behavior by comparison with known inhibitors.

## Background and Rationale

Wound recovery is a complex process involving both the upregulation and downregulation of the same receptors over the course of healing. This creates a challenge for traditional therapies that aim to accelerate wound healing through simple inhibitory or activating activity, since the exact stage of wound repair at the time of treatment cannot be reliably predicted.

However, in the case of **chronic wounds**, the primary issue is fairly well established: the persistent **overactivation** of inflammatory and ECM-degradation receptors. The two key targets here are:

- **Matrix Metalloproteinase-9 (MMP-9)** — involved in extracellular matrix (ECM) degradation
- **Cyclooxygenase-2 (COX-2)** — involved in the inflammatory response

In normal wound healing, both proteins are necessary during the initial stages. In chronic wounds (e.g., diabetic wounds), however, these receptors remain constantly activated, causing excessive ECM breakdown and prolonged over-inflammation, which together impair proper healing.

In contrast, the following receptors are essential for healthy wound healing, particularly in its later stages:

- **Epidermal Growth Factor Receptor (EGFR)** — drives epithelial cell growth
- **Vascular Endothelial Growth Factor Receptor-2 (VEGFR-2)** — drives angiogenesis (new blood vessel formation)

This is the rationale for the project: to identify flavonoid compounds that **selectively inhibit MMP-9 and COX-2** while **avoiding inhibition of EGFR and VEGFR-2**, thereby targeting the pathological drivers of chronic wounds without disrupting the receptors needed for tissue regeneration.

## Objectives

- Screen a select set of 5 flavonoids and their glycones (18 ligands in total), along with known inhibitors and two control inhibitors (**Curcumin** and **Doxycycline**), against all four protein targets
- Compare docked binding affinities across the individual targets
- Evaluate the interaction pattern for each target, using known inhibitors as a reference
- Analyze target selectivity by comparing interaction activity across all four targets

## Protein Targets

| Target | PDB ID | Structure Link |
|---|---|---|
| MMP-9 | 6ESM | [6ESM](https://www.rcsb.org/structure/6ESM) |
| COX-2 | 5IKR | [5IKR](https://www.rcsb.org/structure/5IKR) |
| EGFR | 6VHN | [6VHN](https://www.rcsb.org/structure/6VHN) |
| VEGFR-2 | 3VO3 | [3VO3](https://www.rcsb.org/structure/3VO3) |

All structures were obtained from the [RCSB Protein Data Bank](https://www.rcsb.org/).

## Tools and Software Used

| Tool | Purpose |
|---|---|
| [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) | Environment management |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Ligand structure source |
| [RCSB PDB](https://www.rcsb.org/) | Protein structure source |
| [Open Babel](https://openbabel.org/) | Ligand file format conversion/merging |
| [RDKit](https://www.rdkit.org/) | Hydrogen addition, cheminformatics |
| [Meeko](https://github.com/forlilab/Meeko) | Ligand/receptor PDBQT preparation for AutoDock Vina |
| [PyMOL](https://pymol.org/) | Structure cleaning and visual grid box validation |
| [UCSF Chimera (DockPrep)](https://www.cgl.ucsf.edu/chimera/) | Protein preparation (missing hydrogens, partial charges) |
| [AutoDock Vina](https://vina.scripps.edu/) | Molecular docking (batch mode) |
| [BIOVIA Discovery Studio](https://www.3ds.com/products/biovia/discovery-studio) | Receptor-ligand interaction analysis |

## Workflow

*Simplified Workflow Diagram*
![workflow diagram](./protein%20cleaned/docking%20workflow.png)

### 1. Environment Setup

The project was performed in Ubuntu (WSL). The environment is set up using the provided [`environment.yml`](./environment.yml):

```bash
micromamba create -n <name> -f environment.yml
micromamba activate <name>
```

### 2. Ligand Acquisition

Ligand structures were downloaded from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) in `.sdf` format. For each flavonoid, the most bioactive glycones were identified and selected as ligands.

### 3. Ligand Merging

Individual ligand `.sdf` files were merged per flavonoid using the [Open Babel](https://openbabel.org/) CLI:

```bash
cd <path of ligands folder>
obabel *.sdf -O <flavonoid name>.sdf
```

### 4. Ligand Renaming

Compound names inside each merged `.sdf` file were replaced with common names (for easier downstream processing) using [`sdf_molecule_renamer.py`](./scripts/sdf_molecule_renamer.py):

```bash
python scripts/sdf_molecule_renamer.py
```

### 5. Ligand Preparation

Ligands were prepared using [`ligand_processing.py`](./scripts/ligand_processing.py), which uses **RDKit** to add hydrogen atoms and **Meeko** to convert structures to `.pdbqt` format (assigning atom types and charges):

```bash
python scripts/ligand_processing.py
```

The resulting files are saved to the [`prepared ligands/`](./prepared%20ligands) folder.

### 6. Protein Structure Acquisition

Protein target structures were downloaded from [RCSB PDB](https://www.rcsb.org/). The structures used in this project are listed in the [Protein Targets](#protein-targets) table above.

### 7. Protein Cleaning

Each protein was isolated from solvent molecules and co-crystallized ligands using the **PyMOL** GUI, and the cleaned structures were saved to the [`protein cleaned/`](./protein%20cleaned) folder.

### 8. Reference Ligand Isolation and Processing

The co-crystallized molecules removed in the previous step were retained as reference inhibitor ligands and processed using the same ligand preparation steps described in [Section 5](#5-ligand-preparation).

The receptor and its corresponding co-crystallized ligand were placed together in a folder named after the target's PDB ID. Proteins were then processed with **UCSF Chimera's DockPrep** tool to add missing hydrogen atoms and partial charges.
Example folder: [VEGFR-2 (3VO3)](./3VO3)

### 9. Protein Preparation for Docking

Protein targets were prepared and Vina grid box coordinates generated using Meeko's `mk_prepare_receptor.py`:

```bash
cd <folder of protein>
python $CONDA_PREFIX/bin/mk_prepare_receptor.py --read_pdb <prepared_protein>.pdb -o <protein id> -p -v --box_enveloping <co-crystallized ligand>.sdf --padding 5
```

This command was run for all four protein targets, generating the receptor `.pdbqt` files along with grid box coordinates.

### 10. Grid Box Validation

The prepared protein `.pdb`, isolated co-crystallized ligand, and `.box.pdb` file were loaded together in **PyMOL** to visually confirm that the ligand fits within the grid box, with enough space to accommodate bulkier glycoside ligands.

### 11. Docking Setup

Each protein's `.box.txt` file was updated to include the receptor file path and docking parameters. For each protein, a `ligands` subfolder was created to hold the ligand `.pdbqt` files, and an empty `docked` subfolder was created to store the docking output poses.

### 12. Molecular Docking

Docking was performed using **AutoDock Vina's** batch mode, docking all ligands against each protein target:

```bash
cd <path of individual protein folder>
vina --config <protein id>.box.txt --batch ligands --dir docked 2>&1 | tee docking.log
```

This produces `.pdbqt` files of the docked poses for each ligand in the corresponding `docked/` folder.

### 13. Extracting Binding Affinities

Docking scores were compiled into a `.csv` file using [`extractenergy.py`](./scripts/extractenergy.py):

```bash
python scripts/extractenergy.py
```

### 14. Converting Poses to SDF

Docked poses were converted from `.pdbqt` to `.sdf` format using Meeko's `mk_export.py`. An `sdf files` folder was first created inside each protein's `docked/` directory:

```bash
cd <sdf files folder of a protein>
python $CONDA_PREFIX/bin/mk_export.py ../*.pdbqt --suffix _docked
```

### 15. Pose Selection for Interaction Analysis

Binding affinities were reviewed, and poses within **0.3 kcal/mol** of the top-scoring ligand for each target were selected for interaction analysis. Their `.sdf` files were saved to the [`interaction analysis/`](./interaction%20analysis) folder, under the respective protein subfolder.

### 16. Protein-Ligand Interaction Analysis

Interaction analysis was performed using **BIOVIA Discovery Studio Client**, specifically its Receptor-Ligand Interactions tool:

1. A ligand and protein file were selected, and interactions were visualized as a 2D diagram (with bond lengths enabled in the display settings).
2. Screenshots of the 2D interaction diagrams were saved.
3. This was repeated for **MMP-9** and **COX-2**, and residue interactions were recorded for analysis.
4. Residue interactions were reviewed to identify ligands with predicted inhibitory activity against both MMP-9 and COX-2.
5. These shortlisted ligands were then docked and analyzed against **EGFR** and **VEGFR-2**, and their residue interactions were similarly recorded.

### 17. Selectivity Analysis

Finally, the recorded interactions across all four targets were compared to identify the ligands showing predicted **selective inhibitory activity** — i.e., strong predicted inhibition of MMP-9 and COX-2 without corresponding inhibition of EGFR and VEGFR-2. Detailed results and interaction diagrams for each target are available in the [`interaction analysis/`](./interaction%20analysis) folder.

## Interaction Analysis

This section presents the residue-level interaction analysis for each target, comparing the reference (co-crystallized) inhibitor against the top-ranked docked ligands.

### COX-2 (5IKR)

**Key residues** (see [`interaction analysis/5IKR/observation.txt`](./interaction%20analysis/5IKR/observation.txt) for the full residue-level record):

- Channel entrance residues: **Arg120** and **Tyr355**
- Channel hydrophobic residue: **Val523**
- Catalytic residue: **Ser530**
- **Tyr385** — inhibition here blocks the conversion of arachidonic acid into pro-inflammatory prostaglandins

*Reference ligand bound in the COX-2 active site.*

![reference-5IKR](interaction%20analysis/5IKR/reference.png)

The reference ligand forms a hydrogen bond with the key residue Tyr385, an interaction characteristic of COX-2 inhibitors, along with a π-alkyl interaction with the channel entrance residue Tyr355.

Luteolin, the ligand with the highest binding affinity (**-8.099 kcal/mol**), did not form any of the key interactions described above, despite its favorable docking score.

*Kaempferol bound in the COX-2 active site.*

![kaempferol-5IKR](interaction%20analysis/5IKR/kaempferol.png)

Kaempferol (**-7.811 kcal/mol**) formed two strong hydrogen bonds with Ser530 and a π-cation interaction with Arg120. Together with more than six additional π-alkyl and amide-π interactions contributing to pose stability, this represents a strong overall binding profile. However, kaempferol did not reproduce the characteristic Tyr385 interaction seen in known inhibitors — it is nonetheless still considered a good candidate.

*Isoquercetin bound in the COX-2 active site.*

![isoquercetin-5IKR](interaction%20analysis/5IKR/isoquercetin.png)

Isoquercetin (**-7.495 kcal/mol**) formed two hydrogen bonds with the characteristic residue Tyr385, a π-alkyl interaction with Val523, and two π-cation interactions with Arg120, alongside several further stabilizing contacts. It also showed a weak unfavorable acceptor–acceptor interaction with His90 (3.0 Å); given its low strength, this can be safely disregarded, as it is unlikely to meaningfully disrupt the pose.

### MMP-9 (6ESM)

As MMP-9 is a zinc-dependent metalloproteinase, the key interaction expected of an effective inhibitor is direct metal coordination (chelation) with the catalytic Zn²⁺ ion in the binding pocket — an interaction shown by known inhibitors, including the reference ligand used here.

*Reference ligand chelating the catalytic Zn²⁺ ion in the MMP-9 active site.*

![reference-6ESM](interaction%20analysis/6ESM/reference.png)

None of the docked ligands reproduced this metal-coordination interaction. This is consistent with a known limitation of flavonoids as a chemical class: they lack the functional groups — such as hydroxamic acids, thiolates, or thioesters — required to chelate a zinc atom. As a result, none of the ligands selected for this study are predicted to effectively inhibit MMP-9.

With MMP-9 ruled out, the ligands predicted to inhibit COX-2 were next evaluated against EGFR and VEGFR-2 to assess selectivity.

### EGFR (6VHN)

**Key residues** (see [`interaction analysis/6VHN/observation.txt`](./interaction%20analysis/6VHN/observation.txt) for the full residue-level record):

- Hinge region residue: **Met793**
- Stabilizing salt bridge: **Lys745–Glu762**
- Gatekeeper residue: **Thr790**

*Reference inhibitor bound in the EGFR ATP-binding pocket.*

![reference-6VHN](interaction%20analysis/6VHN/reference.png)

The reference inhibitor forms two strong hydrogen bonds with Met793 and a π-alkyl interaction with Lys745, along with numerous additional stabilizing contacts. These are the key interactions associated with EGFR inhibition — and are therefore interactions this study aims to avoid.

*Kaempferol bound in the EGFR ATP-binding pocket.*

![kaempferol-6VHN](interaction%20analysis/6VHN/kaempferol.png)

Kaempferol did not form any of the key EGFR interactions. It did form strong hydrogen bonds with a few unrelated residues; while not inherently problematic, this is not ideal, as it could still interfere with normal ligand binding within the pocket.

*Isoquercetin bound in the EGFR ATP-binding pocket.*

![isoquercetin-6VHN](interaction%20analysis/6VHN/isoquercetin.png)

Isoquercetin formed a hydrogen bond with Met793, which would ordinarily be a concern for selectivity. However, this pose also carries a critical unfavorable donor–donor clash with the stabilizing residue Lys745 (1.79 Å), making this docking pose difficult to achieve in practice. This makes isoquercetin a potential candidate as well.

### VEGFR-2 (3VO3)

**Key residues** (see [`interaction analysis/3VO3/observation.txt`](./interaction%20analysis/3VO3/observation.txt) for the full residue-level record):

- Kinase domain hinge region: **Cys919**
- DFG motif: **Asp1046, Phe1047, Gly1048**
- Gatekeeper residue: **Val916**

*Reference inhibitor bound in the VEGFR-2 kinase domain.*

![reference-3VO3](interaction%20analysis/3VO3/reference.png)

The reference inhibitor forms two strong hydrogen bonds with Cys919, along with a hydrogen bond and π-cation interaction with Asp1046, though it also carries a relatively weak unfavorable donor–donor clash (2.3 Å). These are the characteristic interactions this study aims to avoid.

*Kaempferol bound in the VEGFR-2 kinase domain.*

![kaempferol-3VO3](interaction%20analysis/3VO3/kaempferol.png)

Kaempferol formed only a single, comparatively weak π-anion interaction with Asp1046, indicating that it could be an ideal candidate.

*Isoquercetin bound in the VEGFR-2 kinase domain.*

![isoquercetin-3VO3](interaction%20analysis/3VO3/isoquercetin.png)

Isoquercetin formed only a single weak π-π stacking interaction with Phe1047, likewise indicating that it could be an ideal candidate.

With the interaction analysis complete, two promising candidates emerged: **kaempferol** and **isoquercetin**.

## Comparative Discussion

It was noted that despite showing a high binding affinity, Luteolin was not an ideal COX-2 inhibitor, as it lacked any of the key residue interactions described above. This demonstrates that docking score alone is not proof of inhibitory activity, and underscores the need for interaction-based analysis of docked poses to correctly identify candidate ligands. This is also why, in standard docking pipelines, all ligands scoring below a practical threshold (typically around -6.5 to -7.0 kcal/mol) are carried forward for interaction analysis, rather than restricting evaluation to only the top-ranked ligand.
At the conclusion of the interaction analysis, two good candidates were identified: kaempferol and isoquercetin.

### Isoquercetin

Isoquercetin — chemically, quercetin 3-*O*-glucopyranoside (also known as quercetin 3-*D*-glucoside) — is a glycosylated derivative of quercetin. It showed a favorable overall predicted inhibition profile against COX-2 owing to its characteristic inhibitor-like interaction with Tyr385. Against EGFR and VEGFR-2, it showed only weak or energetically unfavorable interactions at the key residues, consistent with a favorable selectivity profile. This makes isoquercetin the more promising of the two candidates. As a glycone, it may also gain additional binding stability from its sugar moiety, which can wrap around the protein surface and help anchor the ligand within the binding pocket.

### Kaempferol

Kaempferol showed little favorable interaction with EGFR, which supports its selectivity profile. However, it did not reproduce the characteristic Tyr385 interaction seen in established COX-2 inhibitors. This does not rule out COX-2 inhibition — that could only be confirmed through in vitro assays — but it does mean the predicted activity is less certain than for isoquercetin. Kaempferol is also an aglycone, which may itself be a contributing factor to the docking limitations discussed below (see [Limitations](#limitations)).

## Limitations

- The flavonoids selected for this study were not ideal candidates for inhibiting MMP-9, as they lack the key functional groups required to chelate the catalytic zinc atom.
- COX-2 docking showed a bias toward aglycones over glycones in the binding scores. This is unusual, and is likely a result of not setting flexible side chains at the gatekeeper residues, combined with a grid box too small to accommodate the bulkier glycone ligands without steric hindrance. This is further supported by the fact that most ligands returned only 3–4 docked poses against COX-2, compared to the usual 8–9 poses obtained for the other targets, along with relatively low docking scores overall. This raises the possibility that glycosylated forms of kaempferol — which occur naturally in plants — could in reality show better COX-2 interactions than reflected in these results.

## Proposed Next Steps

- **Improve COX-2 docking performance** by using flexible side-chain docking and a larger grid box, to correct the likely steric bias against glycosylated ligands observed in this study.
- **QSAR-guided optimization** for better binding with MMP-9, or selection of alternative ligand scaffolds better suited to zinc chelation.
- **Molecular dynamics (MD) simulations** to evaluate the conformational stability of the shortlisted docking poses over time, beyond the static snapshot provided by docking alone.
- **MM/PBSA rescoring** of the docked poses following MD, to obtain more reliable, thermodynamically grounded binding free energies than docking scores alone can provide.
- **Enzyme inhibition assays** for the identified candidate ligands, to experimentally validate predicted COX-2 inhibition and target selectivity.
- **Cell-based validation** to confirm functional activity and selectivity in a physiologically relevant system.


