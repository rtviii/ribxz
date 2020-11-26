from typing import Generator, List, Tuple, TypedDict
from Bio.PDB.ResidueDepth import get_surface, ResidueDepth
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.MMCIFParser import FastMMCIFParser
import numpy as np
from asyncio import  run
import argparse
from os import path, stat
import os, sys
from neo4j import GraphDatabase, BoltStatementResult
from Neoget import Neoget
from ResidueMethods import (getLigandsResIds,getAllLigandsForStruct,getResidueNeighbors,getSubchainResidues, ResidueFullIdDict, get_alpha_carbon)
from Ribovision import fetchAlignment
from Structure import (initiateStructureFile)

staticFilesPath = os.getenv( 'STATIC_ROOT' )

async def matchClassToStrand(pdbid:str, banClass:str)->str:
    CYPHER=f"""
    match (r:RibosomeStructure{{rcsb_id: "{pdbid.upper()}"}})-[]-(rp:RibosomalProtein)-[]-(n:NomenclatureClass{{class_id:'{banClass}'}})
    return rp.strand_id"""
    strand_id =Neoget(CYPHER).values()
    if len(strand_id) != 1:
        print('ERROR: Matched none or multiple subchains for a single nomenclature class {} in {}'.format(banClass, pdbid))
        raise NameError
    print(f"Matched strand {strand_id[0][0]} of {pdbid} to nomenclature class {banClass}")
    return strand_id[0][0]




# -----------------------------------------------
def cifparseCli():
    parser   = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='Path to .cif structur file to parse(if not the default defined in the .env)', default='default', nargs='?')
    parser.add_argument('id', type=str,help='PDB-Id of the structure.')
    parser.add_argument('chain', type=str,help='PDB-Id of the subchain.', nargs='*')
    parser.add_argument('-ligsnbr', '--ligandsNeighbors', type=str, dest='ligandsNbr', nargs='?')
    parser.add_argument('--seelig', action='store_true',dest='structLigands')
    parser.add_argument('-lvl', dest='neighborSearchLevel',choices=['C','R'], nargs='?')
    parser.add_argument('-r',   dest='neighborSearchRadius', type=float, nargs='?')
    parser.add_argument('-surf',   dest='surface',  action='store_true')
    args       = parser.parse_args()
    pdbid      = args.id.upper()
    #Download struct if doesn't exist
    if not os.path.exists(os.path.join(staticFilesPath,pdbid,pdbid+".cif")):
        run(initiateStructureFile(pdbid))
    chain         = args.chain
    pathtofile    = args.file
    if args.ligandsNbr:
        level = args.neighborSearchLevel 
        radius = args.neighborSearchRadius
        ligandsNbr  = [str(lig) for lig in args.ligandsNbr.split(',')]
        struct      = fetchStructure(pdbid)
        ligResidues = getLigandsResIds(ligandsNbr, struct)
        for resDict in ligResidues:
            print("\n-------------------------------------------\nLIGAND RESIDUE PROFILE:")
            for item in resDict.__dict__.items():
                print(item)
            nbrsArr = getResidueNeighbors(struct, resDict, radius,level)
            print(f"""\nIn the neighborhood of {radius} Å around ligand [{resDict.chemicalName}]:""")
            for nbr in nbrsArr:
                if level == 'C':
                    nbr:Chain
                    nbr = nbr.get_id()
                    nbr = [nbr, "banClass: {}".format(run( matchStrandToClass(pdbid, nbr) ))] 
                    print(nbr)
                elif level == 'R':
                    print((nbr.get_resname(), nbr.get_id() )," in chain ", nbr.get_parent())
    if args.structLigands:
        for lig in run(getAllLigandsForStruct(pdbid)):
            print(lig)
# -----------------------------------------------

if __name__ == "__main__":
    cifparseCli()