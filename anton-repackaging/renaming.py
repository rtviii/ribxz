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

from ciftools.Neoget import _neoget

prs = FastMMCIFParser(QUIET=True)

io =MMCIFIO()

for pdbid in ['1vy4',]:
    pdbid=pdbid.upper()
    struct:Structure= prs.get_structure(f'{pdbid}', f'{pdbid}.cif')
    for chain in struct[0].child_list:
        strand_id       = chain.id
        nomclass_result = _neoget(f"""match (r:RibosomeStructure{{rcsb_id: "{pdbid.upper()}"}})-[]-(rp:RibosomalProtein)-[]-(n:NomenclatureClass)
        where rp.entity_poly_strand_id  = "{strand_id}" return n.class_id""")
        print(nomclass_result)
        if len(nomclass_result)>0:
            struct[0][strand_id].id = f"{strand_id}[{nomclass_result[0][0]}]"
            print("got nomcalss",nomclass_result[0][0])
        else:
            struct[0][strand_id].id = f"{strand_id}[-]"
            print("got empty")
    
    io.set_structure(struct)
    io.save(f'{pdbid}+.cif')
    


# print(run(matchStrandToClass('3j9m','D')))
# k:Model =struct[0]
# chain:Chain;
# i = 0
# for chain in k.child_list:
#     original = chain.id
#     new      = f"{original}[uL{i}]"
#     i+=1;
#     struct[0][original].id = new

# for cs in struct[0]:
#     print(cs)
# io.set_structure(struct)