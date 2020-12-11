from os import path
from pymol import cmd
from dotenv import load_dotenv
import os,sys


def root(abspath:str, rootname:str)->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    return abspath[:str.find(abspath,rootname) + len(rootname)]




def scoop(pdbid:str, species:int, radius:int):

    envfile= os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    load_dotenv(dotenv_path=envfile)
    STATIC_ROOT  = os.getenv('STATIC_ROOT')
    TUNNELS      = os.getenv('TUNNELS')
    # SCOOP_RADIUS = os.getenv('SCOOP_RADIUS')
    SCOOP_RADIUS = radius

    structfile   = os.path.join(STATIC_ROOT, pdbid, "{}.cif".format(pdbid))
    scoopfile    = os.path.join(TUNNELS,species,pdbid, '{}_{}Ascoop.pdb'.format(pdbid, SCOOP_RADIUS))

    if not os.path.exists(scoopfile):
        os.makedirs(os.path.dirname(scoopfile), exist_ok=True)

    cmd.load(structfile)
    cmd.select('PTC_LIGAND', 'resi 2451')
    cmd.create('TARGET_SELE', '{} w. {} of PTC_LIGAND'.format(pdbid, SCOOP_RADIUS))
    cmd.save(scoopfile, format='pdb',selection='TARGET_SELE')

    print("Saved {}".format(scoopfile))


if __name__ == "__main__":

    pdbid   = sys.argv[1].upper().strip()
    species = sys.argv[2]
    radius  = sys.argv[3]
    scoop(pdbid,species, radius)

