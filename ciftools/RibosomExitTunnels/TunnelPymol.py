#!/usr/bin/env python3

import os,sys,numpy,json
from dotenv import load_dotenv

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


pdibd = sys.argv[1]
struct = LOG._struct(pdbid=pdibd)
struct.render_walls(10)
struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('5j88')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('5u9f')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('5dfe')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('5gag')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('6o9j')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('6s0k')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('6ouo')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')

# struct = LOG._struct('6vwn')
# struct.render_walls(10)
# struct.tunnel_walls.generateReport(write_to_path='.')