import os
from os import error, ftruncate, lseek, mkdir
from typing import List, Tuple
import sys
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.Structure import Structure
from Bio.PDB.mmcifio import MMCIFIO
from Bio import PDB
from asyncio import run
from dotenv import load_dotenv


load_dotenv(dotenv_path='./../.env')
STATIC_ROOT=os.getenv('STATIC_ROOT' )

  
def fetchStructure(pdbid:str, custom_path='default') -> Structure:

    """
    Returns an open PDB.Bio.Structure.Structure object corresponding to <pdbid> from the default repository(specified in the .env)  
    or if custom_path is provided -- from there.
    """

    pathToFile = custom_path if custom_path != 'default' else os.path.join(STATIC_ROOT, pdbid.upper(),  pdbid.upper()+'.cif' )

    if not os.path.exists(pathToFile):
        print(f"File does not exits at the provided path {pathToFile}")
        raise FileNotFoundError(pathToFile) 

    parser:FastMMCIFParser     = FastMMCIFParser(QUIET=True)
    struct:Structure.Structure = parser.get_structure(pdbid.upper(), pathToFile)
    return struct

def fetchChain(pdbid:str, strandid: str, custom_path='default')->Chain:
    pdbid = pdbid.upper()
    struct = fetchStructure(pdbid, custom_path)
    chain:Chain = struct[0][strandid]
    return chain



def splitIntoChains(pdbid:str):
    pdbid = pdbid.upper()
    if len(pdbid) < 2:
        raise IndexError("Provide a valid structure path to parse")
    
    structpath = os.path.join(STATIC_ROOT,pdbid,pdbid+'.cif')
    struct = fetchStructure(pdbid, structpath)[0]

    if not os.path.exists(os.path.join(STATIC_ROOT,pdbid,'CHAINS/')):
        os.mkdir(os.path.join(STATIC_ROOT,pdbid,'CHAINS/'))

    for chain in struct.child_dict.keys():
        destination  = os.path.join(STATIC_ROOT,pdbid,'CHAINS', '{}_STRAND_{}.cif'.format(pdbid,chain) )
        chain = struct[chain]
        io           = MMCIFIO()
        io.set_structure(chain)
        io.save(destination)
        print("Wrote to ",destination)