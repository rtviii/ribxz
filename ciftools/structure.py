from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser
from os import path
import os, sys
from urllib.request import urlretrieve

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

async def initiateStructureFile(pdbid:str)->int:
    staticFilesPath = os.path.join(os.environ['STATIC_ROOT'])
    pdbid = pdbid.upper()
    folder =  os.path.join(staticFilesPath, pdbid)
    structfile = os.path.join(folder, '{}.cif'.format(pdbid))
    
    if not os.path.exists(folder):
        try:
            os.mkdir(folder)
        except: 
            print(f"Failed not create {folder}.")
            raise OSError
    
    urlretrieve('https://files.rcsb.org/download/{}.cif'.format(pdbid), filename=structfile)
    print('Successfully downloaded {}'.format(structfile))
    return 1


