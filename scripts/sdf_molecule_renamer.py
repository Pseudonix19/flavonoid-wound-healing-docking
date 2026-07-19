from rdkit import Chem
import webbrowser
from tkinter import filedialog
import tkinter as tk
import os
root=tk.Tk()
root.withdraw()

inp=filedialog.askopenfilename()
out_dir = filedialog.askdirectory( title="Select Output Folder" )
filename = input("Output filename (without .sdf): ").strip()
output = os.path.join(out_dir, filename + ".sdf")

supplier=Chem.SDMolSupplier(inp)
writer=Chem.SDWriter(output)

for i,mol in enumerate(supplier):
    if mol is None:
        continue

    if mol.HasProp("_Name"):
        curr_name=mol.GetProp("_Name")
    else:
        curr_name="No Name Given"

    if mol.HasProp("PUBCHEM_COMPOUND_CID"):
        cid=mol.GetProp("PUBCHEM_COMPOUND_CID")
    else:
        cid="No CID Given"
        
    print(f"Molecule {i+1}")
    print(f"Current Name: {curr_name}")
    print(f"Current CID: {cid}")
    webbrowser.open_new_tab(f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}")
    
    new_name=input('Enter New Name: ').strip()
    if new_name:
        mol.SetProp("_Name", new_name)

    writer.write(mol)
    
writer.close()
print("Done")
