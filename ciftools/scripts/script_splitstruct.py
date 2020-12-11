import os, sys
from dotenv import load_dotenv
from splitStruct import splitIntoChains

load_dotenv(dotenv_path='./../.env')
STATIC_ROOT=os.getenv('STATIC_ROOT' )

structs = os.listdir(STATIC_ROOT)

for struct in structs:
    if not (os.path.exists(os.path.join(STATIC_ROOT,struct, 'CHAINS'))):
        print(os.path.join(STATIC_ROOT,struct, 'CHAINS')," exists:" , os.path.exists(os.path.join(STATIC_ROOT,struct, 'CHAINS')) )
        splitIntoChains(struct)
    else:
        print(os.path.join(STATIC_ROOT,struct, 'CHAINS') +' exists already.')

