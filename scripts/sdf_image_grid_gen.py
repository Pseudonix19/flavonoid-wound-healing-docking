from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from PIL import ImageDraw
from tkinter import filedialog
import tkinter as tk
import os
root=tk.Tk()
root.withdraw()

inp=filedialog.askopenfilename()
out_dir = filedialog.askdirectory( title="Select Output Folder" )
filename = input("Output filename (without .png): ").strip()
output = os.path.join(out_dir, filename + ".png")

rowsinp=int(input("Enter number of images per row: ").strip())
fntsize=int(input("Enter legend font size: ").strip())

suppl=Chem.SDMolSupplier(inp)

mols=[]
labels=[]
for mol in suppl:
    if mol is None:
        continue
    
    AllChem.Compute2DCoords(mol)
    mols.append(mol)
    
    if mol.HasProp("_Name"):
        name=mol.GetProp("_Name")
    else:
        name="Unknown"
        
    if mol.HasProp("PUBCHEM_COMPOUND_CID"):
        cid=mol.GetProp("PUBCHEM_COMPOUND_CID")
    else:
        cid="Unknown"
    labels.append(f"{name} CID: {cid}")
    
opts=Draw.MolDrawOptions()
opts.legendFontSize=fntsize
opts.padding=0.03
width=int(input("Enter Width of image: ").strip())
height=int(input("Enter Height of image: ").strip())

grid_image=Draw.MolsToGridImage(
    mols,
    molsPerRow=rowsinp,
    subImgSize=(width,height),
    legends=labels,
    drawOptions=opts
    )
draw=ImageDraw.Draw(grid_image)
for x in range(width,grid_image.width,width):
    draw.line((x,0,x,grid_image.height),fill="black",width=2)

for y in range(height,grid_image.height,height):
    draw.line((0,y,grid_image.width,y),fill="black",width=2)
grid_image.save(output)
print("Done")
