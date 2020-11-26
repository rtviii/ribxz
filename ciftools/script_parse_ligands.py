import os, sys
from dotenv import load_dotenv
from ligandDict import parseLigandNeighborhoods

load_dotenv(dotenv_path='./../.env')
STATIC_ROOT=os.getenv('STATIC_ROOT' )



structs = os.listdir(STATIC_ROOT)

for struct in structs:
    print(f"\n-----------------> {struct}")
    struct = struct.upper()
    parseLigandNeighborhoods(struct)
