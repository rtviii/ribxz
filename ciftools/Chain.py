
from asyncio import runners
from ciftools.Structure import fetchStructure
from typing import Generator, List, Tuple, TypedDict
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.MMCIFParser import FastMMCIFParser
import asyncio
import argparse
from os import path, stat
import os, sys
from neo4j import GraphDatabase,BoltStatementResult
from Neoget import Neoget
from ResidueMethods import (getLigandsResIds,getAllLigandsForStruct,getResidueNeighbors,getSubchainResidues, ResidueFullIdDict)
from Ribovision import fetchAlignment
from Structure import (initiateStructureFile)

staticFilesPath = os.path.join(os.getenv['STATIC_ROOT'])

# to get a cif file : curl -O https://files.rcsb.org/download/3j9m.cif 

async def matchClassToStrand(pdbid:str, banClass:str)->str:
    CYPHER=f"""
    match (r:RibosomeStructure{{_PDBId: "{pdbid.upper()}"}})-[]-(rp:RibosomalProtein)-[]-(n:NomenclatureClass{{class_id:'{banClass}'}})
    return rp.strand_id"""
    strand_id =Neoget(CYPHER).values()
    if len(strand_id) != 1:
        print('ERROR: Matched none or multiple subchains for a single nomenclature class {} in {}'.format(banClass, pdbid))
        raise NameError
    print(f"Matched strand {strand_id[0][0]} of {pdbid} to nomenclature class {banClass}")
    return strand_id[0][0]


async def matchStrandToClass(pdbid:str, strand_id:str)->str or None:
    CYPHER=f"""
    match (r:RibosomeStructure{{rcsb_id: "{pdbid.upper()}"}})-[]-(rp:RibosomalProtein)-[]-(n:NomenclatureClass) where rp.strand_id contains "{strand_id}"
    return n.class_id"""
    banClass =Neoget(CYPHER).values()
    if len(banClass) != 1: return None
    else: return banClass


def fetchChain(pdbid:str, strandid: str)->Chain:
    pdbid  = pdbid.upper()
    struct = fetchStructure(pdbid)
    chain:Chain = struct[0][strandid]
    return chain