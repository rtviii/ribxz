#!/usr/bin/env python3
import sys, os, numpy, pandas, json
ciftoolspath= sys.path.append( os.path.dirname(os.path.dirname( os.path.dirname(os.path.realpath(__file__)) ) ))
sys.path.append(ciftoolspath)
from typing import List
from pandas.core.indexes.range import RangeIndex
from dotenv import load_dotenv
# from structure import fetchStructure
# from neoget import Neoget
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import (Chain,Structure,Atom,Residue)





import ciftools