import pandas as pd
import sys, os
ROOT= os.path.dirname(os.path.dirname( os.path.dirname(os.path.realpath(__file__)) ) )
sys.path.append( os.path.dirname(os.path.dirname( os.path.dirname(os.path.realpath(__file__)) ) ))
from pandas.core.indexes.range import RangeIndex
from dotenv import load_dotenv


load_dotenv(dotenv_path=os.path.join(ROOT, '.env'))
AAs          = ["ALA",'ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL']
Nucleotides  = ['A', 'U', 'G', 'C', 'T']

TUNNELS_PATH = os.path.join(ROOT, 'ciftools','TUNNELS')
TUNNEL_LOG   = os.path.join(TUNNELS_PATH, 'TUNNEL_LOG.csv')
frame:pd.DataFrame = pd.read_csv(TUNNEL_LOG)
x      = sys.argv[1]
status = 0
index  = frame[ frame['pdbid']==x ].index.values.astype(int)

if not index.size >0: 
    #doestn exist
    index = len(frame.pdbid) + 1
    frame.loc[index] = [x, status]
else:
    frame.loc[index[0]] = [x, status]
     
frame.to_csv(TUNNEL_LOG,  index=False)