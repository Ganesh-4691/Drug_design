import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen

parser=argparse.ArgumentParser(description='Getting the total number of atoms in smiles as length')
parser.add_argument('--input', type=str, required=True, help='please provide the smiles file(.csv or .dat format)')
parser.add_argument('--data1', type=int, required=True, default='smiles_length.csv', help='output of the smiles_length')
parser.add_argument('--data2', type=int, required=True, help='please provide the smiles file(.csv or .dat format)')
parser.add_argument('--data3', type=int, required=True, default='smiles_length.csv', help='output of the smiles_length')
parser.add_argument('--data4', type=int, required=True, help='please provide the smiles file(.csv or .dat format)')

args=parser.parse_args()

input_file = args.input
o_C1   = args.data1
o_C2   = args.data2
o_C3   = args.data3
o_C4   = args.data4

df = pd.read_csv(input_file)

# define ranges
c1_min, c1_max = o_C1, o_C2
c2_min, c2_max = o_C3, o_C4

f1 = df[
    (df['C1'] >= c1_min) & (df['C1'] < c1_max)
]

f2 = df[
    (df['C2'] >= c2_min) & (df['C2'] < c2_max)
]

filtered = df[
    (df['C1'] >= c1_min) & (df['C1'] < c1_max) &
    (df['C2'] >= c2_min) & (df['C2'] < c2_max)
]

c1 = len(f1)
c2 = len(f2)
count = len(filtered)
p1 = (c1 / len(df)) * 100
p2 = (c2 / len(df)) * 100
percentage = (count / len(df)) * 100

print(f"Molecules in range of C1: {c1}")
print(f"Percentage: {p1:.2f}%")
print(f"Molecules in range of C2: {c2}")
print(f"Percentage: {p2:.2f}%")

print(f"Molecules in range: {count}")
print(f"Percentage: {percentage:.2f}%")

