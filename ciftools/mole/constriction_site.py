from operator import index
import os,sys
import numpy as np
from typing import List
from dotenv import load_dotenv

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')

from ciftools.TunnelLog import (Log)
from ciftools.Neoget import _neoget

# This really shouldn't take more than a few hours.
# I want to be able to retrieve and store coordinates of the constrictions site to initiate mole 
# with them and the PTC as the origin points:
# - retrieve uL22 and uL4 chainnanmes from the db
#
#       - inside the struct, find closest residues ->{ constriction}
#       - get coordinates of the center
#
# - commit each datum to log
# - initiate mole from log

log = Log(os.getenv("TUNNEL_LOG"))

def getConstrictedProteins(pdbid:str)->List[str]:

    cypher="""
    match (n:RibosomeStructure{{rcsb_id:"{pdbid}"}})-[]-(rp:RibosomalProtein) 
    where "uL4" in rp.nomenclature or "uL22" in rp.nomenclature
    return rp.entity_poly_strand_id,rp.nomenclature, n.rcsb_id
    """.format_map({"pdbid":pdbid.upper()})

    response=_neoget(cypher)

    return [*map(lambda x: ( x['rp.nomenclature'][0],x['rp.entity_poly_strand_id'] ), response)]



log.add_column('uL22')
log.add_column('uL4')

for struct in log.all_structs():
    log.update_struct(struct,uL22="AAA",uL4="BBB")


log._write()



#To select rows whose column value equals a scalar, some_value, use ==:
# df.loc[df['favorite_color'] == 'yellow']




# log.add_column('uL22')
# log.add_column('uL4')
# log.add_column('ConstrictionSite')






# print(getConstrictedProteins('5uyl'))







