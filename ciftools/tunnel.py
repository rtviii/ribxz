from argparse import ArgumentParser
import pandas as pd
import os,sys
from dotenv import load_dotenv

def root(abspath:str, rootname:str)->str:
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    return abspath[:str.find(abspath,rootname) + len(rootname)]

ROOT=root(os.path.abspath(__file__), 'ribxz')
sys.path.append(ROOT)
load_dotenv(os.path.join(ROOT,'.env'))

tunenellog = os.getenv('TUNNEL_LOG')


from ciftools.scripts.tunnel_report import TunnelWalls


print(TunnelWalls)

