import os,sys
from pymol import cmd


pdbid=sys.argv[1].upper()

cmd.load(f'./static/{pdbid}/{pdbid}.cif')

cmd.reset()

cmd.ray(300,300)

cmd.png(f"./static/{pdbid}/_ray_{pdbid}.png")



