# ==========================================================================================
# Script: SDF to Descriptor CSV Converter for MPO Scoring
#
# DESCRIPTION:
# This script converts a Maestro-exported SDF file into a CSV file containing molecular
# descriptors necessary for Multi-Parameter Optimization (MPO) scoring.
#
# The following properties are extracted:
# - For RDKit calculations of TPSA:
#     - 'MolWt': Molecular weight
#     - 'MolLogP': logp
#     - 'NumHDonors': Number of hydrogen bond donors
#     - 'NumHAcceptors': Number of hydrogen bond acceptors
#     - 'NumRotatableBonds': Number of rotatable bonds
# - From Maestro/Glide/QikProp SDF tags:
#     - 'docking score': Glide docking score
#     - 'QPlogPo/w': Predicted octanol/water partition coefficient
#     - 'CNS': CNS activity category
#     - 'QPlogBB': Predicted log brain/blood barrier penetration
#
# REQUIRED FORMAT (EXPORTING SDF FROM MAESTRO):
# 1. Go to Table Tool and select the group with the compounds you wish to
#    calculate the MPO score
#       - Make sure you unselect any protein structure and only the 
#       - ligand compounds are selected.
# 2. File → Export Structures:
#     - Format: MDL SD [uncompressed] (.sdf, .sd)
#     - SD Options:
#         - Check "Use display names"
#         - Choose "All" for "Properties to be exported"
#
# USAGE:
# - The script asks the user for:
#     1. The path to the input `.sdf` file
#     2. The name for the output `.csv` file
# - The resulting CSV can then be used in downstream MPO scoring scripts.
#
# OUTPUT:
# - A CSV file saved in the same directory as the input SDF, with one row per compound
#   and all necessary descriptors for MPO scoring.
# ==========================================================================================

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import pandas as pd
import os

# Ask user for SDF file path
sdf_path = input("Enter the full path to your SDF file: ").strip().replace("\\", "/").strip('"')

# Load molecules
suppl = Chem.SDMolSupplier(sdf_path)
mols = [mol for mol in suppl if mol is not None]

# Prepare property extraction
data = []
for mol in mols:
    title = mol.GetProp('_Name') if mol.HasProp('_Name') else "Unknown"

    # Compute descriptors
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)

    # Optional: Add Maestro/QikProp/Glide tags if available
    docking_score = mol.GetProp('docking score') if mol.HasProp('docking score') else ''
    qplogpow = mol.GetProp('QPlogPo/w') if mol.HasProp('QPlogPo/w') else ''
    cns = mol.GetProp('CNS') if mol.HasProp('CNS') else ''
    qplogbb = mol.GetProp('QPlogBB') if mol.HasProp('QPlogBB') else ''
    hbd = mol.GetProp('donorHB') if mol.HasProp('donorHB') else hbd
    hba = mol.GetProp('accptHB') if mol.HasProp('accptHB') else hba
    rot_bonds = mol.GetProp('glide rotatable bonds') if mol.HasProp('glide rotatable bonds') else rot_bonds

    data.append({
        'Title': title,
        'docking score': docking_score,
        'mol MW': mw,
        'donorHB': hbd,
        'accptHB': hba,
        'glide rotatable bonds': rot_bonds,
        'TPSA': tpsa,
        'QPlogPo/w': qplogpow,
        'CNS': cns,
        'QPlogBB': qplogbb
    })

# Create DataFrame
df = pd.DataFrame(data)

# Ask for output CSV name
output_name = input("Enter output CSV filename (without extension): ").strip()
if not output_name:
    output_name = "descriptors_output"

output_dir = os.path.dirname(sdf_path)
output_path = os.path.join(output_dir, f"{output_name}.csv")

df.to_csv(output_path, index=False)
print(f"\nMolecular descriptors saved to: {output_path}")
