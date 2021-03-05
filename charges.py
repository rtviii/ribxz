#!/bin/usr/python3

import sys, os,json
import re
from pathlib import Path
from Bio.PDB.ResidueDepth import residue_depth
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from ciftools.Structure import fetchStructure
from dotenv import load_dotenv
from Bio.PDB import FastMMCIFParser

from functools import  reduce



parser = FastMMCIFParser(QUIET=True)
struct:Structure = parser.get_structure('3j7z', '3J7Z.cif')[0]



with open("3J7Z_TUNNEL_REPORT.json") as infile:
    tunnel = json.load(infile)






strands = tunnel['adjacent_strands']

charge_from_P        =  0
charge_from_GLU_ASP  =  0
charge_from_ARG_LYS  =  0

p_global = []

for strand in strands:
    chain:Chain = struct[strand]
    # Grab all the residues from the chain
    for resobject in strands[strand]:
        if resobject['resname'] in [ 'ASP', 'GLU' ]:
            print("asp/glu")
            charge_from_GLU_ASP-=1;
        if resobject['resname'] in [ 'ARG', 'LYS' ]:
            print("Lys/arg")
            charge_from_ARG_LYS+=1;
    # count the atoms
    residues_in_chain = [* map(lambda x: x['resid'], strands[strand]) ]
    res:Residue
    for res in residues_in_chain:

        # For each residue see the atoms
        residue:Residue  = chain[res]   
        in_atoms = [* residue.get_atoms() ]
        p_atoms = []
        reduce(lambda y,x: y.append(x) if x.get_name()== "P" else "", in_atoms, p_atoms)
        p_global.append(p_atoms)


print(charge_from_GLU_ASP)
print(charge_from_ARG_LYS)
print(len(p_global))









