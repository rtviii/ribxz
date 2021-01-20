from argparse import RawDescriptionHelpFormatter
from copy import Error
from typing import Generator, List, Tuple, TypedDict
from Bio.PDB.Structure import Structure
from nptyping import NDArray
from Bio.PDB.Chain import Chain
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from numpy.core.numeric import full
from Neoget import Neoget


class ResidueFullIdDict(): 
    def __init__(self, res:Residue):
        fullid            = list( res.get_full_id() )
        self.structure    = fullid[0]
        self.model        = fullid[1]
        self.strand_id    = fullid[2]
        self.chemicalName = res.get_resname()
        self.residue_id   = [*fullid[3]][1]




def get_alpha_carbon(res:Residue)->Atom: 
    """Get the first(only) alpha carbon in the residue or an exception"""
    alpha_carbons: List[Atom] = list( filter(lambda internal_atom: internal_atom.get_id() =='CA', res.get_atoms()) )
    if len(alpha_carbons) != 1:
        raise Exception(f"Alpha carbon not found or duplicated in residue {res.get_full_id()}")
    return alpha_carbons[0]

def makeResidueTuple(res:Residue)->Tuple[ str, int, NDArray[3] ]: 
    ac     : Atom       = get_alpha_carbon(res)
    coords : NDArray[3] = ac.get_coord()
    resi   : int        = res.get_id()
    resn   : str        = res.get_resname()
    return tuple(( resn, resi, coords ))

def getSubchainResidues(struct:Structure, subchainId:str) -> Generator[Tuple[ str, int, NDArray[3]], None, None]:

    """
    struct    : the parent struture
    subchainId: pdbx_strand_id of the protein or rna
    ---
    RETURNS a list of tuples for each residue contained in the subchain of the form 
    (the residue's amino acid name, its subchain-local index, the coordinates of its the alpha carbon).
    """

    chain:Chain = struct[0][subchainId]
    residues    = chain.get_residues()

    return (ResidueFullIdDict(res) for res in residues) 

async def getAllLigandsForStruct(pdbid:str)-> List[str]: 

    """
    An asynchronous call to the Neo4j database serving
    ----
    pdbid: the structure of interest
    RETURNS a list of object of the form(chemid, chemname) for each ligand present in the struct
    """

    pdbid=pdbid.upper()
    CYPHER_STRING=f"""
    match (r:RibosomeStructure {{`rcsb_id`:"{pdbid}"}})-[]-(l:Ligand)
    return collect(l.chemicalId)
    """

    result = Neoget(CYPHER_STRING)
    return result.values()[0]

def getLigandsResIds(ligChemIds:List[str], struct: Structure)->List[ResidueFullIdDict]:
    """Returns a list of dictionaries specifying each ligand as a residue inside a given struct."""
    ligandResidues: List[Residue] = list(filter(lambda x: x.get_resname() in ligChemIds, list( struct.get_residues() )))
    return [ ResidueFullIdDict(res) for res in ligandResidues ]
    
def getLigandResIds(ligchemid:str, struct: Structure, dict_or_res:bool='dict')->List[ResidueFullIdDict]:
    """Returns a list of dictionaries specifying each ligand as a residue inside a given struct."""
    ligandResidues: List[Residue] = list(filter(lambda x: x.get_resname() == ligchemid, list( struct.get_residues() )))
    if dict_or_res == 'dict':
        return [ ResidueFullIdDict(res) for res in ligandResidues ]
    else:
        return ligandResidues

def getResidueNeighbors(struct:Structure, resid: ResidueFullIdDict, radius:float, level:str = 'R', all_levels=False)-> List[Residue] or List[Chain] or List[Chain or Residue]:

    """
    struct: opened cif structure
    resid: A dictionary containing the residue's identifiers (see associated class)
    radiue: radius of neighborhood to check
    level: Atom, Residue, Chain : one of A/R/C
    """

    if level.upper() not in ['A', 'R', 'C']:
        print('Level has to be one of A R C')
        raise Error

    parentStrand      :Chain         = struct[resid.model][resid.strand_id]
    residuesOfInterest:List[Residue] = list(filter(lambda x : x.get_full_id()[3][1] == resid.residue_id, parentStrand.child_list))

    ligandAtoms = residuesOfInterest[0].get_atoms()

    coords      = list( map(lambda atom: atom.get_coord(), ligandAtoms) )
    ns          = NeighborSearch(list( struct.get_atoms() ))

    yield ns.search(coords[0],radius,level.upper())
     

def getLigandNeighbors (struct: Structure, constituents: List[Residue], rad :float):
    ligatoms = [  r.get_atoms() for r in constituents ]
    print(ligatoms)



def getResidueNeighborsMinimal(struct:Structure, resid: ResidueFullIdDict, radius:float, level:str = 'R', all_levels=False)-> List[Residue] or List[Chain] or List[Chain or Residue]:

    """
    struct: opened cif structure
    resid: A dictionary containing the residue's identifiers (see associated class)
    radiue: radius of neighborhood to check
    level: Atom, Residue, Chain : one of A/R/C
    """

    if level.upper() not in ['A', 'R', 'C']:
        print('Level has to be one of A R C')
        raise Error

    parentStrand      :Chain         = struct[resid.model][resid.strand_id]
    residuesOfInterest:List[Residue] = list(filter(lambda x : x.get_full_id()[3][1] == resid.residue_id, parentStrand.child_list))

    ligandAtoms = residuesOfInterest[0].get_atoms()

    coords      = list( map(lambda atom: atom.get_coord(), ligandAtoms) )
    ns          = NeighborSearch(list( struct.get_atoms() ))

    yield ns.search(coords[0],radius,level.upper())
     