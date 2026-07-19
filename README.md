# Structure-Based Analysis of Flavonoid Selectivity for Chronic Wound Healing Targets Using Molecular Docking and Interaction Analysis

*In simpler words: analyzing flavonoids to find which ones could potentially help in chronic wound healing.*
*(Written to be accessible to student peers who want to follow or reproduce this pipeline.)*

---

## Table of Contents

- [Abstract](#abstract)
- [Background and Rationale](#background-and-rationale)
- [Objectives](#objectives)
- [Protein Targets](#protein-targets)
- [Repository Structure](#repository-structure)
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

Ligands were prepared using [`ligand_processing.py`](./scripts), which uses **RDKit** to add hydrogen atoms and **Meeko** to convert structures to `.pdbqt` format (assigning atom types and charges):

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

Docking scores were compiled into a `.csv` file using [`extractenergy.py`](./scripts):

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

---

*This README documents the full pipeline from ligand/protein preparation through docking and selectivity analysis, and is intended to help peers reproduce or build on this workflow.*
