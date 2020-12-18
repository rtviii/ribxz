

# Khanh's implementation
# xyz=cmd.get_coords(’Tunnel5’,1)
# r=[]
# cmd.iterate_state(1,’Tunnel’,’r.append(vdw)’,space=locals(),atomic=0)
# python
# from pymol import stored
# np.savetxt(’tunnel_coordinates.txt’,xyz,fmt=’\%.2f’)
# np.savetxt(’tunnel_radius.txt’,r,fmt=’\%.2f’)
# python end

import os,sys
from dotenv import load_dotenv

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')

from ciftools.TunnelLog import Log

from ciftools.Structure import fetchStructure
from ciftools.scripts.splitStruct import fetchChain
from typing import List
from Bio.PDB.Residue import Residue
from pymol import cmd


RADIUS      = os.getenv('SCOOP_RADIUS')
TUNNELS     = os.getenv('TUNNELS')
STATIC_ROOT = os.getenv('STATIC_ROOT')

def get_ptc_residues(struct)->List[Residue]:

    def belongs_to_ptc(x:Residue):
        return str(x.get_id()[1]) in ["2055","2056","2451","2452","2507","2506"]

    PTC_residues = filter(belongs_to_ptc, [*struct.get_residues()]) 

    return [* PTC_residues ]

pdbid         = sys.argv[1].upper()
log           = Log(os.getenv('TUNNEL_LOG'))
struct        = fetchStructure(pdbid)
record        = log.get_struct(pdbid)
species       = str( record.taxid.values[0] )

structpath    = os.path.join(STATIC_ROOT,pdbid,"{}.cif".format(pdbid))
scoopsavepath = os.path.join(TUNNELS,species,pdbid,'{}_{}Ascoop.pdb'.format(pdbid,RADIUS))
if not os.path.exists(os.path.dirname( scoopsavepath )):
    os.makedirs(os.path.dirname( scoopsavepath ))


cmd.load(structpath)
x = cmd.select(f'resi 2506')
cmd.select(f'br. {pdbid} w. {RADIUS} of \'sele\'')
cmd.save( scoopsavepath, 'sele' )
print("Saved to {}".format(scoopsavepath))






# if __name__ =='__main__':

#     print(f"Executing {__file__} as a standalone script")

#     parser = argparse.ArgumentParser(f"Argparser for {__file__}", add_help=""" p3 mole/scoop.py -radius 50 -p ./MOLEtrials/4ug0/pymol/4ug0.cif -c L5  -r 4452 4UG0""")
#     parser.add_argument('-radius','-scoop_radius', dest='radius',help="the radius of matter to extract around a given residue")

#     # parser.add_argument('-o','--origin',dest='origin', help='origin of the sphere to scoop matter from.', required=False)
#     parser.add_argument('-p,','--cifpath', help='Relative path to the cif file of the molecule to extract matter from.')
#     parser.add_argument('-c','--chain')
#     parser.add_argument('-r','--residue')
#     parser.add_argument('pdbid', help='pdbid')
#     args    = parser.parse_args()
#     path    = args.cifpath
#     pdbid   = args.pdbid.lower()

#     # origin  = args.origin
#     radius  = args.radius
#     chain   = args.chain
#     residue = args.residue



# def scoop():
#     print("got args \n",args)
#     cmd.load(path)
#     x = cmd.select(f'c. {chain} and resi {residue}')
#     print(f'br. {pdbid} w. {radius} of \'sele\'')
#     cmd.select(f'br. {pdbid} w. {radius} of \'sele\'')
#     cmd.save( os.path.join(os.path.dirname(path), f'scoop_{pdbid}_{radius}.pdb'), 'sele' )


# scoop()