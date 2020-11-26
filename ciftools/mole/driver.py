#!/usr/bin/env python
import json
import os
import xml
import numpy
from pymol import cmd
import pymol
from makeconfig import make_input_config
from cli_parser import makeparser
import asyncio
from dotenv import load_dotenv
import numpy as np

load_dotenv(dotenv_path='./.env')


# sys.path.append(os.path.join(sys.path[0], 'mole2'))
# sys.path.append(sys.path[0])
# path_moleexe = os.path.join(sys.path[0], 'mole2', 'mole2.exe')

# # path_moleconfig = os.path.join(sys.path[0],'input.xml')
# path_moleconfig = os.path.join(sys.path[0], 'input-multiorigin.xml')
# path_pdbstruct = os.path.join(os.path.dirname(sys.path[0]), 'static', 'pdb-structs', 'radiusobject.pdb')



if __name__ == "__main__":

    """
    So here is how this works:

    pass MOLE's arguments as sysargs to this script or predefine them below. 
    (An additional set of defaults is already written in makeconfig.)
    An xml file is then created with parameters for a particular run of mole.
    Don't forget to specify the input STRUCTURE file and OUTPUT were mole would write.
    """
    # load cli options
    moleparser    =  makeparser()

    load_dotenv(dotenv_path='./../../.env')

    MOLE_EXE       =  os.getenv('MOLE_EXECUTABLE')
    ECOLI_TUNNELS  =  os.getenv("ECOLI_TUNNELS")

    args       = moleparser.parse_args()
    #get rid of the unutilized arg options

    args        =  filter((lambda kvpair: None not in kvpair), vars(args).items())
    args        =  dict(args)

    pdbid       =  args['PDBID'].upper()
    if not os.path.exists(os.path.join(ECOLI_TUNNELS,pdbid)):
        os.makedirs(os.path.join(ECOLI_TUNNELS,pdbid) )

    inputconfigpath  =  os.path.join(ECOLI_TUNNELS, pdbid, '{}_moleinput.xml'.format(pdbid) )
    inputstructpath  =  os.path.join(ECOLI_TUNNELS, '{}_SCOOP60A_RESI2451.pdb'.format(pdbid))
    outpath          =  os.path.join(ECOLI_TUNNELS, pdbid)

    args['Input']         =  inputstructpath
    args['Output']        =  outpath
    args['ConfigPath']    =  inputconfigpath

    pymol.cmd.load(inputstructpath)
    pymol.cmd.select('PTC', 'resi 2506')
    cords = pymol.cmd.get_coords('PTC', 1)

    pts = []
    for cord in cords:
        pt =   numpy.around(cord, 0)
        print(pt)
        pts.append([ "{},{},{}".format(pt[0], pt[1], pt[2]) ])



    args['Points']              =  pts

    args['ProbeRadius']         =  "10"
    args['InteriorThreshold']   =  "0.8"
    args['BottleneckRadius']    =  "1"
    args['SurfaceCoverRadius']  =  "10"
    args['OriginRadius']        =  "5"

    args['exports'] = "t"

    #Sensisble inputs 
    asyncio.run(make_input_config(args))
    os.system("mono {} {}".format(MOLE_EXE, inputconfigpath))
