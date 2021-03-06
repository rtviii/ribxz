from typing import List, Tuple
import sys,os
from os import path
import pandas as pd
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
import json
from asyncio import run
from dotenv import load_dotenv


def fetchStructure(pdbid:str, custom_path='default') -> Structure:
    """
    Returns an open PDB.Bio.Structure.Structure object corresponding to <pdbid> from the default repository(specified in the .env)  
    or if custom_path is provided -- from there.
    
    """
    pathToFile = custom_path if custom_path != 'default' else path.join(os.getenv('STATIC_ROOT'), pdbid.upper(),  pdbid.upper()+'.cif' )

    if not path.exists(pathToFile):
        print(f"File does not exits at the provided path {pathToFile}")
        raise FileNotFoundError(pathToFile) 
    parser:FastMMCIFParser     = FastMMCIFParser(QUIET=True)
    struct:Structure.Structure = parser.get_structure(pdbid.upper(), pathToFile)
    return struct

# uri         = os.getenv('NEO4J_URI')
# authglobal  = (os.getenv('NEO4J_USER'),os.getenv('NEO4J_PASSWORD'))
# current_db  = os.getenv('NEO4J_CURRENTDB')

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    root=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(root)
    load_dotenv(os.path.join(root,'.env'))


root_self('ribxz')
STATIC_ROOT = os.getenv('STATIC_ROOT')
from ciftools.Neoget import _neoget

class ResidueFullIdDict(): 
    def __init__(self, res:Residue):
        fullid            = list( res.get_full_id() )
        self.structure    = fullid[0]
        self.model        = fullid[1]
        self.strand_id    = fullid[2]
        self.chemicalName = res.get_resname()
        self.residue_id   = [*fullid[3]][1]

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

def getLigandResIds(ligchemid:str, struct: Structure, dict_or_res:bool='dict')->List[ResidueFullIdDict]:
    """Returns a list of dictionaries specifying each ligand as a residue inside a given struct."""
    ligandResidues: List[Residue] = list(filter(lambda x: x.get_resname() == ligchemid, list( struct.get_residues() )))
    if dict_or_res == 'dict':
        return [ ResidueFullIdDict(res) for res in ligandResidues ]
    else:
        return ligandResidues

def filterIons(entry):
    if "ion" in entry['name'].lower():
       print("Filtered ION:", entry['name'])
       return False
    else:
        return entry['id']

async def matchStrandToClass(pdbid:str, strand_id:str)->str:
    CYPHER="""match (r:RibosomeStructure{{rcsb_id: "{}"}})-[]-(rp:RibosomalProtein{{entity_poly_strand_id:"{}"}})-[]-(n:NomenclatureClass)
    return n.class_id""".format(pdbid.upper(), strand_id)

    # CYPHER="""match (r:RibosomeStructure{rcsb_id:"4U3N"})-[]-(rp:RibosomalProtein)-[]-(n:NomenclatureClass) where rp.entity_poly_strand_id contains "s8"
    # return n""".format()

    resp = _neoget(CYPHER)
    print(CYPHER)
    print("RESP", resp)

    if len(resp) > 0:
        return resp[0]
    else:
        return None

class ResidueId():

    def __init__(self, res:Residue):
        fid             = list(res.get_full_id())
        self.resn    = fid[3][0]
        self.struct     = fid[0]
        self.strand_id  = fid[2]
        self.residue_id = [*fid[3]][1]

    def __eq__(self, other):
        return self.residue_id == other.residue_id and self.strand_id == other.strand_id

    def __hash__(self):
        return hash(( 'strand_id',self.strand_id,'residue_id', self.residue_id ))

    def asdict(self):
        return {
            "resn"     : self.resn,
            "strand_id": self.strand_id,
            "resid"    : self.residue_id,
            "struct"   : self.struct}

def addBanClass(x:ResidueId):
    """Tag a ResidueId with a BanClass"""

    profile = x.asdict()
    bc      = run(matchStrandToClass(profile[ 'struct' ],profile[ 'strand_id' ]))
    print(f"Residue {profile['resn']}/{profile['resid']}", " belongs to strand",profile['strand_id'], ". Identified as Ban class",  bc)
    profile['banClass'] = bc
    return profile

def getLigandNbrs(resids: List[Residue], struct:Structure):
    ns   = NeighborSearch(list( struct.get_atoms() ))
    nbrs = []
    for r in resids:
        # a ligand consists of residues
        resatom = r.child_list[0]
        #  each residue has an atom plucked at random
        for nbr in ns.search(resatom.get_coord(), 5,level='R'):
            # we grab all residues in radius around that atom and extend the list with those
            nbrs.extend([* nbr])
    filtered = [] 
    for neighor in nbrs:
        present = 0
        for constit in resids:
            if ResidueId( constit ) == ResidueId( neighor ):
                present = 1
        if present == 0:
            filtered.append(ResidueId(neighor))
    return [ * map(lambda x: addBanClass(x) ,  set(filtered) ) ]

def parseLigandNeighborhoods(pdbid:str):
    """
    @pdibd is the 4-letter rcsb id.
    """
    pdbid=pdbid.upper()

    entry:List     = _neoget("""match (l:Ligand)-[]-(r:RibosomeStructure{{rcsb_id:"{pdbid}"}}) 
    return {{struct: r.rcsb_id, ligs: collect({{ id:l.chemicalId, name: l.chemicalName }})}}""".format_map({ "pdbid":pdbid }))[0]

    if len(entry)  == 0:
        print(f"No ligands found for {pdbid} the DB. Exiting..")
        return
    else:
        print("Received ligands for {}: ".format(pdbid), entry)

    ligandsResponse:List[str] = entry[0]['ligs']
    ligandIds = [* map(lambda x : filterIons(x), ligandsResponse) ]
    ligandIds = [i for i in ligandIds if i] 
    
    pathtostruct = os.path.join(STATIC_ROOT,pdbid,'{}.cif'.format(pdbid))

    """Iterate over"""
    for x in ligandIds:
        savepath = os.path.join(STATIC_ROOT, pdbid, 'LIGAND_{}.json'.format(x))

        #! These structures and ligand either take forever to render or fail silently. Why?

        if x in ["A"]:
            continue

        # if pdbid in ['4U3N'] or x in ['A', 'OHX', ]:
        #     print("Skipping problematic {}".format(pdbid))
        #     continue
        # if pdbid in ['5TGM'] or x in ['A', 'OHX', 'LEU']:
        #     continue

        # if os.path.exists(savepath):
        #     print(savepath, " already exists. Skipping rendering.")
        #     continue
        # else:

        struct              = fetchStructure(pdbid, pathtostruct)
        print("Parsing residues of {}".format(x))
        asresiudes = getLigandResIds(x, struct, 'res')
        internals  = [addBanClass( ResidueId(residue) ) for residue in asresiudes ]
        nbrs       = getLigandNbrs(asresiudes, struct)

        for nbr in nbrs:
            run(matchStrandToClass(nbr[ 'struct' ],nbr[ 'strand_id' ]))
        ligprofile = {'constituents': internals,'nbrs':         nbrs}

        with open(savepath, 'w') as json_file:
            json.dump(ligprofile,json_file)
            print(f'Wrote to {savepath}')



if __name__ == "__main__":
    pdbid = sys.argv[1]
    parseLigandNeighborhoods(pdbid=pdbid)
