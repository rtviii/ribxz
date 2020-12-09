#!/usr/bin/env python3

import os,sys,numpy,json
from dotenv import load_dotenv
import matplotlib.pyplot  as plt

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')

from ciftools.RibosomExitTunnels.TunnelLog import ( TunnelRecord, TunnelWalls, Log)

LOG = Log(os.getenv('TUNNEL_LOG'))

struct = LOG._struct('5j7l')
struct.render_walls()
struct.tunnel_walls.generateReport(write_to_path='.')

struct = LOG._struct('5j88')
struct.render_walls()
struct.tunnel_walls.generateReport(write_to_path='.')

