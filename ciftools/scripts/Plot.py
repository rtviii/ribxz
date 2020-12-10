import os,sys,json
import numpy as np
from typing import Iterator, List
from numpy.core import numeric
from numpy.core.fromnumeric import _ptp_dispatcher
import pandas as pd
from Bio.PDB.Residue import Residue
from dotenv import load_dotenv
import matplotlib.pyplot  as plt
import matplotlib.patches  as patches
from mpl_toolkits import mplot3d
from numpy.core.records import record

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    ROOT = os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

root_self('ribxz')

from  ciftools.scripts.TunnelLog import (Log,TunnelRecord,TunnelWalls, get_CA_or)

PDBID                  = sys.argv[1].upper()
LOG                    = Log(os.getenv('TUNNEL_LOG'))
STRUCT_REPORT_PATH     = os.path.join('.', f'{PDBID}_TUNNEL_REPORT.json')
ProteinsColorgenerator = iter([ 'red','purple','cyan','orange','green',"gray",'yellow','magenta','pink' ])

report = {}
with open(STRUCT_REPORT_PATH) as infile:
    report=  json.load(infile)

record      = LOG._struct(PDBID)
record.plot_radius()
ptcresidues = record.tunnel_walls.get_ptc_residues()
tunneldf    = record.get_total_df()

def locate_res_in_df(rescoords:List[float],tunnel_df:pd.DataFrame): 
    respos   = np.array(rescoords)
    rowid    = -1
    currdist = 10000000

    for index,row in tunnel_df.iterrows():

        centerlinepos = np.array( ( row['X'],row['Y'], row['Z'] ) )
        dist = np.linalg.norm(centerlinepos-respos)
        if dist <  currdist:
            currdist = dist
            rowid    = index
    
    return tunnel_df.loc[rowid]


proteinsLegend=[]
for prot in report['proteins'].keys():
    curr_color   = next(ProteinsColorgenerator)
    banName = report['nomMap'][prot][0] if len( report['nomMap'][prot] ) >0 else prot
    proteinsLegend.append(patches.Patch(color=curr_color,label=banName))

    chain_coordinates   = [*map(lambda res: res['rescoord'],report['proteins'][prot])]
    chain_corresp_radii = [*map(lambda respos: locate_res_in_df(respos, tunneldf), chain_coordinates)]
    chain_corresp_radii = [*map(lambda row: ( row['Distance'],row['Radius'] ),chain_corresp_radii )]
    dist                = [x[0]for x in chain_corresp_radii]
    radii               = [x[1]for x in chain_corresp_radii]
    plt.scatter(dist,radii,c=curr_color)


ptc_res_locs      = [* map(lambda res: get_CA_or(res).get_coord().tolist(),ptcresidues) ]
ptc_corresp_radii = [*map(lambda respos: locate_res_in_df(respos, tunneldf), ptc_res_locs)]
ptc_corresp_radii = [*map(lambda row: ( row['Distance'],row['Radius'] ),ptc_corresp_radii )]
dist              = [x[0]for x in ptc_corresp_radii]
radii             = [x[1]for x in ptc_corresp_radii]
plt.scatter(dist,radii, c='blue')


axs=plt.axes()
axs.legend(handles=[ *proteinsLegend, patches.Patch(color='blue',label='PTC') ])
axs.set_title(PDBID)
plt.show()
    

