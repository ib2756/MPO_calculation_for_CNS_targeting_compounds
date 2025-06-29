# CNS MPO Scoring and Partial Agonist Screening Toolkit
This repository contains a four-stage Python-based pipeline for evaluating ligand candidates targeting central nervous system (CNS) proteins—specifically focused on the discovery of dopamine D3 receptor (DRD3) partial agonists. The workflow integrates cheminformatics tools to systematically extract molecular descriptors, 
apply CNS-specific multi-parameter optimization (MPO) scoring, and perform benchmark-based filtering of candidate molecules across multiple receptor conformations.
The scoring strategy implemented in 1_CNS_MPO_calculator.py was developed with blood-brain barrier (BBB) permeability and CNS drug-likeness in mind. 
The selected molecular properties (e.g., QPlogPo/w, QPlogBB, TPSA, HBD/HBA, rotatable bonds) and their respective desirability ranges were based on pharmacokinetic profiles of known CNS-penetrant drugs. The flat-top Gaussian desirability function ensures a smooth, penalty-aware evaluation of each parameter, 
while normalizing docking scores provides fair comparison of binding strength across molecules. While this MPO scoring is essential for early prioritization, it is insufficient on its own for identifying partial agonists—which must exhibit optimal interactions in both the active and inactive conformations of the GPCR target. To address this, the pipeline includes 2_Partial_Agonist_matchmaker.py and 3_Partial_Agonist_Postprocessor.py. 

---

## Scripts Overview

### `0_SDF_to_CSV_converter.py`
**Purpose**: Converts Maestro-exported `.sdf` structure files into `.csv` files containing essential molecular descriptors (e.g., TPSA, MW, QPlogPo/w, QPlogBB, CNS).

- Automatically computes TPSA using RDKit
- Parses Glide and QikProp properties from SD tags
- Output CSV is ready for MPO scoring

---

###  `1_CNS_MPO_calculator.py`
**Purpose**: Calculates Multi-Parameter Optimization (MPO) scores using flat-top Gaussian desirability functions for CNS drug-like profiling.

- Normalizes Glide docking scores
- Supports histogram visualization of score distributions
  - MPO < 0.4 compounds are highlighted
  - Docking score spread is annotated with best/worst values

---

### `2_Partial_Agonist_matchmaker.py`
**Purpose**: Cross-matches compounds between two receptor states (e.g., active vs. inactive DRD3), keeping only those outperforming a benchmark (e.g., cariprazine) in both.

This script compares two MPO-scored ligand datasets—typically derived from docking against active and inactive conformations of the same GPCR (e.g., DRD3). 
It retains only compounds that outperform a benchmark ligand in both states. By doing so, it isolates ligands with balanced receptor engagement, 
a critical requirement for partial agonist behavior.

- Input: Two MPO-scored CSVs from different receptor conformations
- Filters based on minimum benchmark MPO
- Output: Matched compounds with MPO above threshold in both files

---

### `3_Partial_Agonist_Postprocessor.py`
**Purpose**: Computes average/delta MPO and normalized docking scores for compound pairs and filters those surpassing a benchmark.

This step ensures that selected candidates are not only better in MPO score than the reference drug in both binding states, but also show minimal disparity between 
the conformations—consistent with partial agonist binding profiles that stabilize intermediate receptor states.

- Sorts by average MPO and average docking separately
- Filters compounds that outperform the benchmark in **both** metrics
- Appends benchmark compound rows for comparison

---

## Example Outputs

- `Sorted_by_Avg_MPO.csv`: Ranked compounds by average MPO
- `Sorted_by_Avg_normDocking.csv`: Ranked compounds by average normalized docking score
- `Above_cariprazine_Compounds.csv`: Final candidates outperforming benchmark on both metrics

---

## Requirements

- Python 3.x
- RDKit
- pandas
- matplotlib (for histogram visualization)

Install RDKit via Conda (to enable full features of the RDKit):
(Use 'conda install -c rdkit rdkit' on the terminal to install instead of 'pip install rdkit')

---

## Recommended Workflow

#### A. Docking & QikProp (in Maestro)
Perform Glide docking and QikProp property prediction on your ligand set using PDB structures of both the active and inactive conformations of the target GPCR protein.

#### B. Export Structures from Maestro
Export the ligand structures as .sdf files from Maestro’s Table tool.
*** Be sure to: Enable "Use display names" and Select "All properties" for export

#### C. Generate Descriptors
Run *0_SDF_to_CSV_converter.py* to convert the exported .sdf file into a .csv file containing all required molecular descriptors.

#### D. MPO Scoring
Use *1_CNS_MPO_calculator.py* to compute MPO scores based on CNS drug-likeness criteria and normalize the docking scores.

#### E. Match Active/Inactive States
Run *2_Partial_Agonist_matchmaker.py* to align compounds across active and inactive conformations, filter out those that underperform a benchmark (e.g., cariprazine), and retain only those with superior dual-state performance.

#### F. Postprocess and Rank Candidates
Execute *3_Partial_Agonist_Postprocessor.py* to calculate average and delta MPO/docking scores, rank candidates, and output top-performing partial agonist leads.

---

## Author
In Young Bae of New York City College of Technology developed this toolkit as part of a CNS-targeted partial agonist drug design project, with the ChatGPT 4o model providing substantial support during the coding implementation stage.
