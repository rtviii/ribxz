#!/usr/bin/env python3

import asyncio
import os,sys
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

    from ciftools.TunnelScripts.TunnelLog import Log

    """
    So here is how this works:

    pass MOLE's arguments as sysargs to this script or predefine them below. 
    (An additional set of defaults is already written in makeconfig.)
    An xml file is then created with parameters for a particular run of mole.
    Don't forget to specify the input STRUCTURE file and OUTPUT were mole would write.

    """
    # load cli options
    moleparser = makeparser()

    load_dotenv(dotenv_path='/home/rxz/dev/ribxz/.env')

    MOLE_EXE     = os.getenv('MOLE_EXECUTABLE')
    TUNNELS      = os.getenv("TUNNELS")
    SCOOP_RADIUS = os.getenv("SCOOP_RADIUS")

    args         = moleparser.parse_args()
    args        =  filter((lambda kvpair: None not in kvpair), vars(args).items())
    args        =  dict(args)

    pdbid       = args['PDBID'].upper()
    load_dotenv('/home/rxz/dev/ribxz/.env')

    log     = Log(os.getenv("TUNNEL_LOG"))
    record  = log.get_struct(pdbid)

    species         = str( int( record.taxid.values[0] ) )
    constriction    = record.constriction_coord.values[0]

    inputconfigpath = os.path.join(TUNNELS, species, pdbid, '{}_moleinput.xml'.format(pdbid) )
    inputstructpath = os.path.join(TUNNELS, species, pdbid, '{}_{}Ascoop.pdb'.format(pdbid, SCOOP_RADIUS))
    outpath         = os.path.join(TUNNELS, species, pdbid)

    args['Input']         =  inputstructpath
    args['Output']        =  outpath
    args['ConfigPath']    =  inputconfigpath

    # pymol.cmd.load(inputstructpath)
    # pymol.cmd.select('PTC', 'resi 2506')
    # cords = pymol.cmd.get_coords('PTC', 1)

    origins = []

    # exits   = []

    # for ptccord in cords:
    #     origins.append([ "{},{},{}".format(ptccord[0], ptccord[1], ptccord[2]) ])

    cord = log.get_struct(pdbid)['PTC_average'].values[0]
    origins.append([cord])

# ---------------------------------
    #original
    # args['Points']             = origins
    # args['ProbeRadius']        = "10"
    # args['InteriorThreshold']  = "0.8"
    # args['BottleneckRadius']   = "1"
    # args['SurfaceCoverRadius'] = "10"
    # args['OriginRadius']       = "5"

    #After Khanh's feedback, tunnel cutoffs on the edge of the scoop.
    args['Points']             = origins
    args['ProbeRadius']        = "12"
    args['InteriorThreshold']  = "1.25"
    args['SurfaceCoverRadius'] = "10"
    args['OriginRadius']       = "5"
    args['BottleneckRadius']   = "1"
    args['exports']            = "t"
# ---------------------------------
    asyncio.run(make_input_config(args))
    os.system("mono {} {}".format(MOLE_EXE, inputconfigpath))