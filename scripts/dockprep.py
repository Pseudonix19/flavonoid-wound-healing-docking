from rdkit import Chem
from rdkit.Chem import AllChem
from tkinter import filedialog
from meeko import MoleculePreparation, PDBQTWriterLegacy
preparator=MoleculePreparation()
import tkinter as tk
root=tk.Tk()
root.withdraw()
from pathlib import Path

def get_mols(selbatch=0,seldir=0):
	if selbatch:
		inp=Path(filedialog.askdirectory( title="Select SDFs directory"))
		(inp.parent / "Result").mkdir(parents=True,exist_ok=True)
	else:
		inp=Path(filedialog.askopenfilename())
	if seldir:
		out_dir = Path(filedialog.askdirectory( title="Select Output Folder" ))
	else:
		out_dir=inp.parent / inp.stem if not selbatch else inp / "Result"
		out_dir.mkdir(parents=True, exist_ok=True)
	if selbatch:
		sdfs=list(inp.glob("*.sdf"))
	else:
		sdfs=[inp]
	return sdfs,out_dir

def prep_mol(mol,i):
	mol=Chem.AddHs(mol, addCoords=True)
	cid = mol.GetProp("PUBCHEM_COMPOUND_CID") if mol.HasProp("PUBCHEM_COMPOUND_CID") else "unknown_cid"
	setup=preparator.prepare(mol)[0]
	pdbqt_data=PDBQTWriterLegacy.write_string(setup)
	nameprop=mol.GetProp("_Name").split()[0] if mol.HasProp("_Name") else f"unknown_name_{i}"
	filename=f"{nameprop}_{cid}"
	return pdbqt_data, filename, nameprop
