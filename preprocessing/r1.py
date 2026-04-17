from rdkit import Chem

sup = Chem.SDMolSupplier(
    "D_new.sdf",
    removeHs=False,
    strictParsing=False
)

writer = Chem.SDWriter("D_new1.sdf")

for mol in sup:
    if mol:
        try:
            Chem.SanitizeMol(mol)
            writer.write(Chem.RemoveHs(mol))
        except:
            pass

writer.close()

print("Finished")
