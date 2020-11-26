from pymol import cmd
from dotenv import load_dotenv
import os,sys



###
ecoli_structures = ["7K00", "5AFI", "3J9Y", "3J9Z", "3JA1", "3JCJ", "4UY8", "5GAD", "5GAE", "5GAG", "5GAH", "5JTE", "5JU8", "5L3P", "5LZA", "5LZD", "5LZE", "5MDV", "5MDW", "5MDZ", "5MGP", "5U9F", "5U9G", "5WDT", "5WE4", "5WE6", "5WF0", "5WFS", "6B4V", "6BOH", "6BOK", "6C4I", "6DNC", "6ENF", "6ENJ", "6ENU", "6GC0", "6GWT", "6GXM", "6GXN", "6GXO", "6HRM", "6I0Y", "6I7V", "6O9J", "6OGF", "6OGI", "6ORE", "6ORL", "6OSK", "6OSQ", "6OT3", "6OTR", "6OUO", "6OXA", "6OXI", "6PJ6", "6Q95", "6Q97", "6Q9A", "6S0K", "6U48", "6VU3", "6VWL", "6VWM", "6VWN", "6VYQ", "6VYR", "6VYS", "6WD0", "6WD1", "6WD2", "6WD3", "6WD4", "6WD5", "6WD7", "6WD8", "6WD9", "6WDA", "6WDD", "6WDE", "6WDK", "6WDM", "6WNT", "6WNV", "6WNW", "6X7F", "6X7K", "6XDQ", "6YSR", "6YSS", "6YST", "6YSU"]
###

load_dotenv(dotenv_path='./../.env')
STATIC_ROOT    =  os.getenv('STATIC_ROOT')
ECOLI_TUNNELS  =  os.getenv('ECOLI_TUNNELS')





for pdbid in ecoli_structures:

    pdbid = pdbid.upper()
    structfile  =  os.path.join(STATIC_ROOT, pdbid, "{}.cif".format(pdbid))
    scoopfile   =  os.path.join(ECOLI_TUNNELS, '{}_SCOOP60A_RESI2451.pdb'.format(pdbid))


    # PYmol

    cmd.load(structfile)
    cmd.select('PTC_LIGAND', 'resi 2451')
    cmd.create('TARGET_SELE', '{} w. 60 of PTC_LIGAND'.format(pdbid))
    cmd.save(scoopfile, format='pdb',selection='TARGET_SELE')
    print("Saved {}".format(scoopfile))