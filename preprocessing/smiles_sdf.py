import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def csv_smiles_to_sdf(csv_filepath, smiles_column_name, output_sdf_path):
    print(f"Reading data from {csv_filepath}...")
    
    # 1. Read the CSV file
    try:
        df = pd.read_csv(csv_filepath)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return
    
    # 2. Check if the SMILES column actually exists
    if smiles_column_name not in df.columns:
        print(f"Error: Could not find a column named '{smiles_column_name}' in your CSV.")
        print(f"Available columns are: {list(df.columns)}")
        return

    # 3. Set up the single SDF writer
    writer = Chem.SDWriter(output_sdf_path)
    
    success_count = 0
    fail_count = 0

    # 4. Loop through the SMILES row by row
    for index, row in df.iterrows():
        smi = str(row[smiles_column_name]).strip()
        
        # Skip empty rows
        if not smi or smi == 'nan':
            continue
            
        mol = Chem.MolFromSmiles(smi)
        
        if mol is not None:
            try:
                # Add Hydrogens and generate 3D coordinates
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                
                # --- BONUS ---
                # Save all the other columns from your CSV as properties in the SDF!
                for col in df.columns:
                    mol.SetProp(str(col), str(row[col]))
                
                # Write to the single file
                writer.write(mol)
                success_count += 1
            except Exception as e:
                fail_count += 1
        else:
            fail_count += 1
            
        # Print progress so you know it's working
        if (index + 1) % 500 == 0:
            print(f"Processed {index + 1} rows... ({success_count} succeeded, {fail_count} failed)", end="\r")

    # Close the file when done
    writer.close()
    print(f"\n✅ Finished! Successfully converted {success_count} molecules.")
    print(f"❌ Failed: {fail_count}")
    print(f"📁 All molecules saved to one file: {output_sdf_path}")

if __name__ == "__main__":
    
    # --- CHANGE THESE THREE VARIABLES BEFORE RUNNING ---
    
    CSV_FILE = "data.csv"        # The name of your CSV file
    SMILES_COLUMN = "SMILES"              # The exact header name of the column containing the SMILES strings
    OUTPUT_SDF = "300.sdf"      # The name of the single output SDF file
    
    # -------------------------------------------------

    csv_smiles_to_sdf(CSV_FILE, SMILES_COLUMN, OUTPUT_SDF)
