import dockprep
from rdkit import Chem

selbatch=input("Batch prep multiple sdfs? Leave blank for no ")
seldir=input("Select output directory? Leave blank for autogeneration ")
sdfs,out_dir=dockprep.get_mols(selbatch,seldir)

total=0
failed=0
completed=0

for sdf in sdfs:
	supplier=Chem.SDMolSupplier(sdf)

	for i,mol in enumerate(supplier):
		total+=1
		if mol is None:
			failed+=1
			continue

		data,filename,molname=dockprep.prep_mol(mol,i)

		if not data[1]:
			print(f"Failed {molname}")
			print(data[2])
			failed+=1
			continue

		with open(out_dir / f"{filename}.pdbqt", "w") as f:
			f.write(data[0])
		completed+=1

print(f"Processed: {total}")
print(f"Succeeded: {completed}")
print(f"Failed: {failed}")
