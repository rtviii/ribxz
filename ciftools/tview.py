from pymol import cmd
import os, sys
from dotenv import load_dotenv


load_dotenv(dotenv_path='./../.env')



def tview(pdbid:str):
    pdbid          =  pdbid.upper()

    MOLE_EXE       =  os.getenv('MOLE_EXECUTABLE')
    ECOLI_TUNNELS  =  os.getenv("ECOLI_TUNNELS")
    TUNNEL_SCRIPT  =  os.path.join(ECOLI_TUNNELS, pdbid, 'pymol', 'complex.py')

    inputstructpath  =  os.path.join(ECOLI_TUNNELS, '{}_SCOOP60A_RESI2451.pdb'.format(pdbid))
    cmd.load(inputstructpath)
    cmd.color('gray', 'all')
    cmd.select('PTC', 'resi 2055 or resi 2056 or resi 2451 or resi 2452 or resi 2507 or resi 2506')
    cmd.color('blue', 'PTC')
    cmd.run(TUNNEL_SCRIPT)

    cmd.hide('everything', 'Tunnels')
    cmd.show('mesh', 'Tunnels')
    cmd.reset()



def twrite(pdbid, args):
    pdbid= pdbid.upper()
    args           =  args.split(' ')
    ECOLI_TUNNELS  =  os.getenv("ECOLI_TUNNELS")
    CHOICE_TXT     =  os.path.join(ECOLI_TUNNELS, pdbid, 'csv', 'choices.txt')

    f = open(CHOICE_TXT, 'w')
    f.write(",".join(args))
    f.close()

cmd.extend("tview", tview)
cmd.extend("twrite", twrite)