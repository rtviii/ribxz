import sys
from typing import List
import os 
import numpy
from pandas.core.indexes.range import RangeIndex
from structmethods import fetchStructure
from dotenv import load_dotenv
import pandas as pd
from Neoget import Neoget
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import (Chain,Structure,Atom,Residue)
import json


envp = './../.env'
load_dotenv(dotenv_path=envp)


AAs         = ["ALA",'ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL']
Nucleotides = ['A', 'U', 'G', 'C', 'T']
STATIC_ROOT    =  os.getenv('STATIC_ROOT')
ECOLI_TUNNELS  =  os.getenv("ECOLI_TUNNELS")

if __name__=="__main__":

    pdbid='';
    tunneln=[]

    if sys.argv[1] == '--ecoli':
        print("Working with ecoli presets.")
        pdbid = sys.argv[2].upper()
        print("Got pdbid :", pdbid)
        TUNNELS_PATH  =  os.path.join(ECOLI_TUNNELS,pdbid,'csv')
        CHOICE_TXT    =  os.path.join(TUNNELS_PATH, 'choices.txt')
        f = open(CHOICE_TXT)
        for csv_number in f.read().split(","):
            tunneln.append(csv_number)

    else:
        pdbid   = sys.argv[1].upper()
        tunneln = sys.argv[2].split(',')

        TUNNELS_PATH   =  os.path.join(STATIC_ROOT,pdbid,'TUNNEL','csv')

    if len(pdbid) != 4 or len(tunneln) < 1:
        print("Enter a valid RCSB Id followed by comma-sepated tunnel ids: [4UG0 1,2]")
        exit(1)
else:
    print("Please run me as a module. [ -m ]")
    exit(1)



tunnelfiles = []
for t in tunneln:
    tpath = os.path.join(TUNNELS_PATH, f'tunnel_{t}.csv')
    tunnelfiles.append(tpath)

total = pd.DataFrame()
for tunnelpath in tunnelfiles:
    inst  = pd.read_csv(tunnelpath)
    xyzr  = inst[["FreeRadius", "X","Y","Z"]]
    total = total.append(xyzr, ignore_index=True)

struct = fetchStructure(pdbid)
atoms  = list(struct.get_atoms())
ns     = NeighborSearch(atoms,bucket_size=10)

class TunnelWalls:
    def __init__(self, pdbid:str) -> None:
        self.struct = pdbid.upper()
        self.rna                = {}
        self.rps                = {}
        self.nomenclatureMap    = {}
        self.other              = []
        self.adjacentRnaStrands = []
        self.adjacentRPStrands  = []
        self.radius             = []
        self.ligands            = []
        self.rescount           = 0

    def getProteinResidues(self):
        return self.rps

    def getRnaResidues(self):
        return self.rna

    def addResidue(self, res:Residue.Residue)->None:
        parentStrand = res.get_parent().get_id()
        if res.get_resname() not in [*AAs, *Nucleotides] and res not in self.ligands:
            self.ligands.append(res)
        if parentStrand not in self.adjacentRnaStrands and parentStrand not in self.adjacentRPStrands:
            response = Neoget("""
            match (n{{entity_poly_strand_id:"{parentStrand}"}})-[]-\
            (r:RibosomeStructure{{rcsb_id:"{struct}"}}) \
            return {{type: n.entity_poly_polymer_type, nomenclature:n.nomenclature}};""".format_map({
                "parentStrand":parentStrand,
                "struct":self.struct
            }))
            profile  = response.values()[0][0]

            if profile['type'] == 'RNA':
                self.adjacentRnaStrands.append(parentStrand)
                self.rna[parentStrand] = []

            if profile['type'] == 'Protein':
                self.adjacentRPStrands.append(parentStrand)
                self.rps            [parentStrand] = []
                self.nomenclatureMap[parentStrand] = profile['nomenclature']

        if parentStrand in self.adjacentRnaStrands:
            if res in self.rna[parentStrand] :
                None
            else:
                self.rna[parentStrand].append(res)
        elif parentStrand in self.adjacentRPStrands:
            if res in self.rps[parentStrand]:
                None
            else:
                self.rps[parentStrand].append(res)

        self.rescount +=1

    def consumeMoleDataframe(self,df:pd.DataFrame, radius:float):
        """
        Takes a dataframe which must contain X,Y,Z columns assuming the centerline of the tunnel,
        although other columns can be present. Aimed at MOLE's mergd csv results.
        Iterates over each row and applies neighbor search on each focus, appends non-redundant residues
        to appropriate registry on the object. 
        """
        self.radius = radius
        def getVicinity(row):
            x = row['X'];
            y = row['Y'];
            z = row['Z']
            res:List[Residue.Residue] = ns.search(numpy.array([x,y,z]), radius,level='R')
            for nbr in res:
                self.addResidue(nbr)
        df.apply(getVicinity, axis=1)
    
    def generateReport(self, pathtowrite):
        def getResInfo(res:Residue.Residue, nomMap, polytype:str ): 
            parent                                                : Chain.Chain     = res.get_parent()
            parentStrand                                          : str             = parent.get_id()
            resid                                                 : int             = res.get_id()[1]
            resname                                               : str             = res.get_resname()
            nom   = nomMap[parentStrand] if polytype == "Protein" else None 
            islig = False if resname.upper() in [*Nucleotides, *AAs] else True
            return {"strand": parentStrand, "resid" : resid, "resname": resname, "polytype": polytype, "nom": nom, "isligand": islig}
        
        protwall  =  {}
        rnawall   =  {}

        presentLigands = [getResInfo(l,nomMap=self.nomenclatureMap, polytype="Other") for l in tws.ligands]

        for tpl in self.getProteinResidues().items():
            protwall[ tpl[0] ] = [getResInfo(x,nomMap=self.nomenclatureMap, polytype="Protein") for x in tpl[1]]
        for tpl in self.getRnaResidues().items():
            rnawall[ tpl[0] ] = [getResInfo(x,nomMap=self.nomenclatureMap, polytype="RNA") for x in tpl[1]]

    
        rprt           = {
            "pdbid"        :  pdbid,
            "probeRadius"  :  self.radius,
            "rna"          :  rnawall,
            "proteins"     :  protwall,
            "ligands"      :  presentLigands,
            "nomMap"       :  self.nomenclatureMap

        }

        with open(pathtowrite, 'w') as writable:
            json.dump(rprt,writable)


    
        

    
tws = TunnelWalls(pdbid=pdbid)
tws.consumeMoleDataframe(total, 10)
reportPath = os.path.join(STATIC_ROOT,pdbid,'TUNNEL',f"{pdbid}_TUNNEL_REPORT.json")
if not os.path.exists(os.path.dirname( reportPath )):
    os.makedirs(os.path.dirname( reportPath ))
    
report =tws.generateReport(reportPath)





