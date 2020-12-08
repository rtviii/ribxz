import os,sys,numpy,json
from dotenv import load_dotenv
import matplotlib.pyplot  as plt
from mpl_toolkits import mplot3d

def root_self(rootname:str='')->str:

    """Returns the rootpath for the project if it's unique in the current folder tree."""
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))
root_self('ribxz')

from  ciftools.RibosomExitTunnels.TunnelLog import (Log,TunnelRecord,TunnelWalls)

PDBID='3JCD'
LOG         = Log(os.getenv('TUNNEL_LOG'))
STRUCT_REPORT_PATH = os.path.join(os.getenv('STATIC_ROOT'),PDBID,'TUNNEL', f'{PDBID}_TUNNEL_REPORT.json')


report = {}
with open(STRUCT_REPORT_PATH) as infile:
    report=  json.load(infile)

record_3jcd = LOG._struct('3jcd')
df          = record_3jcd.get_total_df()
xs          = df['X'].tolist()
ys          = df['Y'].tolist()
zs          = df['Z'].tolist()


fig = plt.figure()
ax  = plt.axes(projection='3d')
ax.scatter3D(xs,ys,zs)
plt.show()



