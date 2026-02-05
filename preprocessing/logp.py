import argparse
from rdkit import Chem
from rdkit.Chem import Crippen
import pandas as pd

# --- Step 1: Read SMILES from file ---
parser = argparse.ArgumentParser(description="Compute logP values from SMILES file")
parser.add_argument("--input", type=str, required=True, help="Path to your data file (.dat or .csv)")
parser.add_argument("--output", type=str, default="output_with_logp.csv", help="Output CSV file name")
args = parser.parse_args()

file_path = args.data
output_path = args.output

with open(file_path, "r") as f:
    smiles_list = [line.strip() for line in f if line.strip()]

# --- Step 2: Compute logP for each molecule ---
data = []
for smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        logp = Crippen.MolLogP(mol)
        data.append((smi, logp))
    else:
        data.append((smi, None))

# --- Step 3: Store in DataFrame ---
df = pd.DataFrame(data, columns=["SMILES", "logP"])
print(df.head())

# --- Step 4: Save results if you want ---
df.to_csv("smiles_with_logp.csv", index=False)
print("✅ Saved file: smiles_with_logp.csv")

