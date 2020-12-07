from pymol import cmd
import os, sys
from dotenv import load_dotenv



species      = '83333'
def root(abspath:str, rootname:str)->str:
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    return abspath[:str.find(abspath,rootname) + len(rootname)]


def tview(pdbid:str):

    load_dotenv(dotenv_path="/home/rtviii/dev/ribxz/.env")

    pdbid        = pdbid.upper()
    TUNNELS      = os.getenv("TUNNELS")
    SCOOP_RADIUS = os.getenv("SCOOP_RADIUS")
    TUNNEL_SCRIPT   = os.path.join(TUNNELS,species, pdbid, 'pymol', 'complex.py')
    inputstructpath = os.path.join(TUNNELS,species, pdbid, '{}_{}Ascoop.pdb'.format(pdbid,SCOOP_RADIUS))


    print("input", inputstructpath)

    cmd.load(inputstructpath)
    cmd.color('gray', 'all')
    cmd.select('PTC', 'resi 2055 or resi 2056 or resi 2451 or resi 2452 or resi 2507 or resi 2506')
    cmd.color('blue', 'PTC')
    if not os.path.exists(TUNNEL_SCRIPT):
        print("TUNNEL SCRIPT NOT FOUND.")
        twrite(pdbid, "0")
        return
    cmd.run(TUNNEL_SCRIPT)

    cmd.hide('everything', 'Tunnels')
    cmd.show('mesh', 'Tunnels')
    cmd.reset()



def twrite(pdbid, args):
    # Anticipating the missing tunnels, usually would write PDBID, [1 2 3] to choose tunnel.
    # In the case of faulty struct/absence of tunnels : write PDBID, 0

    # Write -1 for a chain present inside the tunnel

    pdbid         = pdbid.upper()
    args          = args.split(' ')

    TUNNELS      = os.getenv("TUNNELS")
    CHOICE_TXT    = os.path.join(TUNNELS,species, pdbid,  'tunnels-results.txt')

    # if len(args) == 1 and int( args[0] ) == 0:
    # then it failed thats what's up. just write it

    f = open(CHOICE_TXT, 'w')
    f.write(",".join(args))
    f.close()
    print(f"Wrote results to {CHOICE_TXT}")
    cmd.delete('all')
    print("Cleared all, pymol.")

cmd.extend("tview", tview)
cmd.extend("twrite", twrite)