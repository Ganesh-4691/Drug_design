import argparse
import pandas as pd
from rdkit import Chem

parser=argparse.ArgumentParser(description='Getting the total number of atoms in smiles as length')
parser.add_argument('--input', type=str, required=True, help='please provide the smiles file')
parser.add_argument('--output', type=str, default='smiles_length.csv', help='output of the smiles_length')
args=parser.parse_args()

input_file = args.data
output_file   = args.output

with open(input_file, "r") as f:
    smiles_list=[line.strip() for line in f if line.strip()]
results=[]
for smi in smiles_list:
    mol=Chem.MolFromSmiles(smi)
    if mol:
        atom_count=mol.GetNumAtoms()
        results.append((smi,atom_count))
    else:
        results.append((smi, None))


df=pd.DataFrame(results, columns=["smiles", "length"])
print(df.head)

df.to_csv(output_file,index=False)
print("length of the smiles are generated")
