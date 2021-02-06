from Bio.PDB.Entity import Entity
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.mmcifio import MMCIFIO
from dotenv import load_dotenv
import os, sys
from asyncio import run

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    root=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(root)
    load_dotenv(os.path.join(root,'.env'))

root_self('ribxz')



from ciftools.Ligand import parseLigandNeighborhoods
