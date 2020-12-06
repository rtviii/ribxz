#!/usr/bin/env python

import asyncio

import pymol
import json,os,xml,numpy,sys
sys.path.append(os.path.dirname( os.path.dirname( os.path.dirname(os.path.realpath(__file__)) ) ))
from ciftools.mole.cli_parser import makeparser
from ciftools.mole.makeconfig import make_input_config
from dotenv import load_dotenv



def root(abspath:str, rootname:str)->str:
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    return abspath[:str.find(abspath,rootname) + len(rootname)]

if __name__ == "__main__":

    """
    So here is how this works:

    pass MOLE's arguments as sysargs to this script or predefine them below. 
    (An additional set of defaults is already written in makeconfig.)
    An xml file is then created with parameters for a particular run of mole.
    Don't forget to specify the input STRUCTURE file and OUTPUT were mole would write.

    """
    # load cli options
    moleparser = makeparser()
    envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')

    load_dotenv(dotenv_path=envfile)

    MOLE_EXE     = os.getenv('MOLE_EXECUTABLE')
    TUNNELS      = os.getenv("TUNNELS")
    SCOOP_RADIUS = os.getenv("SCOOP_RADIUS")
    args         = moleparser.parse_args()
    #get rid of the unutilized arg options

    args        =  filter((lambda kvpair: None not in kvpair), vars(args).items())
    args        =  dict(args)

    pdbid   = args['PDBID'].upper()
    species = args['taxid']
    
    inputconfigpath = os.path.join(TUNNELS, species, pdbid, '{}_moleinput.xml'.format(pdbid) )
    inputstructpath = os.path.join(TUNNELS, species, pdbid, '{}_{}Ascoop.pdb'.format(pdbid, SCOOP_RADIUS))
    outpath         = os.path.join(TUNNELS, species, pdbid)

    args['Input']         =  inputstructpath
    args['Output']        =  outpath
    args['ConfigPath']    =  inputconfigpath

    pymol.cmd.load(inputstructpath)
    pymol.cmd.select('PTC', 'resi 2506')
    cords = pymol.cmd.get_coords('PTC', 1)

    pts = []
    for cord in cords:
        pt =   numpy.around(cord, 0)
        pts.append([ "{},{},{}".format(pt[0], pt[1], pt[2]) ])

    args['Points']             = pts
    args['ProbeRadius']        = "10"
    args['InteriorThreshold']  = "0.8"
    args['BottleneckRadius']   = "1"
    args['SurfaceCoverRadius'] = "10"
    args['OriginRadius']       = "5"

    args['exports'] = "t"

    #Sensisble inputs 
    asyncio.run(make_input_config(args))
    possible = [os.path.join(TUNNELS, species, pdbid, '{}_{}Ascoop.pdb'.format(pdbid, x)) for x in ['80','60','70']]

    for scoop in possible:
        args['Input']         =  scoop
        try:
            os.system("mono {} {}".format(MOLE_EXE, inputconfigpath))
            break
        except:
            continue
