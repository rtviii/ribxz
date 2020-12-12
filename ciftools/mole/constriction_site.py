import os,sys
from dotenv import load_dotenv

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')

from ciftools.TunnelLog import (Log)


