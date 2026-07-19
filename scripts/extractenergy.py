from pathlib import Path
from tkinter import filedialog
import tkinter as tk
root=tk.Tk()
root.withdraw()
import pandas as pd

inp=Path(filedialog.askdirectory( title="Select PDBQTs directory"))
out=Path(filedialog.askdirectory(title="Select Results.csv directory"))
pdbqts=list(inp.glob("*.pdbqt"))
sep=float(input("Energy difference between best energy and least useful energy: "))

results=[]
for file in pdbqts:
	nameinfo=file.stem.split("_")
	if nameinfo[1].isdigit():
		with open(file) as f:
			lines=f.readlines()
		val=[]
		for line in lines:
			if "REMARK VINA RESULT:" in line:
				val.append(line.split()[3])
		val=list(map(float, val))
		best=min(val)
		pose=val.index(best) +1
		record = {
			"Name": nameinfo[0],
			"CID": nameinfo[1],
			"Affinity": best,
			"Pose": pose
			}
		results.append(record)
df=pd.DataFrame(results)
df=df.sort_values("Affinity").reset_index(drop=True)
print(df)
print()
fileprefix=out.parent.name
df.to_csv(out / f"{fileprefix}_Results.csv", index=False)
max=min(df["Affinity"])
top=df[df["Affinity"] <= max + sep]
top=top[["Name", "Affinity"]]
print(top)
top.to_csv(out / "sdf files" / f"{fileprefix}_Useful.csv", index=False)
