#!/usr/bin/env python3

import os,sys,json,numpy
from dotenv import load_dotenv
def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))
root_self('ribxz')

from typing import List
from dotenv import load_dotenv
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import (Chain,Structure,Atom,Residue)
import pandas as pd
from ciftools.structure import fetchStructure
from ciftools.neoget import _neoget


AAs          = ["ALA",'ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL']
Nucleotides  = ['A', 'U', 'G', 'C', 'T']

class TunnelRecord:

    def __init__(self,data_path:str, pdbid:str, taxid:int, moleoutputs:List[int], comments:List[str]=[""]):
        self.__struct_data_loc    =  data_path
        self.pdbid                =  pdbid
        self.taxid                =  taxid
        self.csv_choices          =  moleoutputs
        self.comments             =  comments
        self.total_df             =  pd.DataFrame()
        self.computed_properties  =  {"cumulative_distance":  0}
        """cumulative_distance | proteins_map | """
        self.concat_mole_csvs()

    def concat_mole_csvs(self):

        total             =  pd.DataFrame()
        distance          =  0
        mole_chosen_csvs  =  [os.path.join(self.__struct_data_loc,'csv', 'tunnel_{}.csv'.format(x)) for x in self.csv_choices]

        for f in mole_chosen_csvs:
            tunnel_instance   =  pd.read_csv(f)
            xyzr              =  tunnel_instance[['Distance','FreeRadius', 'X','Y','Z']]
            rows              =  len(xyzr)
            distance         +=  xyzr.iloc[rows-1]['Distance']
            total.append(xyzr)
            
        self.computed_properties['cumulative_distance'] = distance
        self.total_df =total

class TunnelWalls:
    def __init__(self, pdbid:str, structure: Structure) -> None  :  
        self.structure  =  structure
        self.pdbid      =  pdbid.upper()

        self.rna                 =  {}
        self.rps                 =  {}

        self.nomenclatureMap     =  {}
        self.other               =  []
        self.adjacentRnaStrands  =  []
        self.adjacentRPStrands   =  []
        self.radius              =  []
        self.ligands             =  []

        self.rescount            =  0


    def getProteinResidues(self):
        return self.rps

    def getRnaResidues(self):
        return self.rna

    def addResidue(self, res:Residue.Residue)->None:
        parentStrand = res.get_parent().get_id()
        if res.get_resname() not in [*AAs, *Nucleotides] and res not in self.ligands:
            self.ligands.append(res)

        if parentStrand not in self.adjacentRnaStrands and parentStrand not in self.adjacentRPStrands:

            response = _neoget("""
            match (n{{entity_poly_strand_id:"{parentStrand}"}})-[]-(r:RibosomeStructure{{rcsb_id:"{struct}"}}) \
            return {{type: n.entity_poly_polymer_type, nomenclature:n.nomenclature}};""".format_map({
                "parentStrand": parentStrand,
                "struct"      : self.struct
            }))
            try:
                profile = response[0]
                print(profile)
                if profile['type'] == 'RNA':
                    self.adjacentRnaStrands.append(parentStrand)
                    self.rna[parentStrand] = []

                if profile['type'] == 'Protein':
                    self.adjacentRPStrands.append(parentStrand)
                    self.rps            [parentStrand] = []
                    self.nomenclatureMap[parentStrand] = profile['nomenclature']
            except:
                pass

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
        print(self.rescount)

    def consumeMoleDataframe(self,df:pd.DataFrame, radius:float):

        """
        Takes a dataframe which must contain X,Y,Z columns assuming the centerline of the tunnel,
        although other columns can be present. Aimed at MOLE's mergd csv results.
        Iterates over each row and applies neighbor search on each focus, appends non-redundant residues
        to appropriate registry on the object. 
        """
        atoms        =  list(self.structure.get_atoms())
        ns           =  NeighborSearch(atoms,bucket_size=10)
        self.radius  =  radius

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

        report           = {
            "pdbid"        :  self.struct,
            "probeRadius"  :  self.radius,
            "rna"          :  rnawall,
            "proteins"     :  protwall,
            "ligands"      :  presentLigands,
            "nomMap"       :  self.nomenclatureMap

        }

        with open(pathtowrite, 'w') as writable:
            json.dump(report,writable)

class Log:

    """
    Log file itself is the interface to the csv files produced by MOLE.
    """

    def __init__(self, path:str)->None  :
        """Logging utility for keeping track of the ribosomal tunnels,
        cosuming and concatenating MOLE csv outputs,
        miscellaneous comments. \
        Initiate from a ${PROJECT_ROOT}/ciftools/TUNNELS/TUNNEL_LOG.csv
        _______________________________________________________________
        An assumption is made that the TUNNELS_LOG.csv resides in the 
        top-level directory which contains all the other tunnel files.
        """
        self.log:pd.DataFrame   = pd.read_csv(path)
        self.__tunnels_path     = os.path.dirname(path)

    def _all_structs(self):
        yield self.log['pdbid'].tolist()
        
    def _struct(self,pdbid:str)->TunnelRecord:
        pdbid  =  pdbid.upper()
        row    =  self.log.loc[self.log['PDBID'] ==pdbid]
        return TunnelRecord(
        os.path.join(self.__tunnels_path,str(row['TaxId'][0]),row['PDBID'][0]),
        row['PDBID'][0],
        row['TaxId'][0],
        row['MoleStatus'][0].split(','))

log            =  Log(os.getenv('TUNNEL_LOG'))
struct_record  =  log._struct('3jcd');

TunnelWalls('3jcd', fetchStructure('3jcd'))





