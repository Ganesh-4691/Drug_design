import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen

parser=argparse.ArgumentParser(description='Getting the total number of atoms in smiles as length')
parser.add_argument('--input', type=str, required=True, help='please provide the smiles file(.csv or .dat format)')
parser.add_argument('--output', type=str, required=True, default='smiles_length.csv', help='output of the smiles_length')
args=parser.parse_args()

input_file = args.input
output_file   = args.output

with open(input_file, "r") as f:
    smiles_list=[line.strip() for line in f if line.strip()]
results=[]
for smi in smiles_list:
    mol=Chem.MolFromSmiles(smi)
    if mol:
        atom_count=mol.GetNumAtoms()
        logp_value=Crippen.MolLogP(mol)
        mr_value=Crippen.MolMR(mol)
        results.append((smi,atom_count, logp_value, mr_value))
    else:
        results.append((smi, None))


df=pd.DataFrame(results, columns=["smiles", "atom_count", "logP", "MR"])
print(df.head)

df.to_csv(output_file,index=False)
print("length of the smiles are generated")

df_original = pd.read_csv(output_file)
atom_num = "atom_count"
df_filtered = df_original[df_original[atom_num] <=16]
df_filtered.to_csv("filtered_smiles_max16atoms.csv", index=False)

print("Done")
