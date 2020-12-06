from argparse import ArgumentParser
from typing import List
import pandas as pd
import os,sys
from dotenv import load_dotenv
from ciftools.scripts.tunnel_report import TunnelWalls

def root_self(rootname:str='')->str:
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

root_self('ribxz')
STATIC_ROOT = os.getenv('STATIC_ROOT')
TUNNEL_LOG  = os.getenv('TUNNEL_LOG')
TUNNELS     = os.getenv('TUNNELS')



AAs         = ["ALA",'ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL']
Nucleotides = ['A', 'U', 'G', 'C', 'T']

#



def load_tunnels( pdbid:str, species:int , tunnel_ids:List[int]=[] ):
    """Load tunnel_ids as recorded in the log into the dataframe"""
    MOLE_CSVS =os.path.join(TUNNELS, species, pdbid,'csv')



    with open(TUNNEL_RESULTS) as chosentunnels:
        chosentunnels = chosentunnels.read()

        df:pd.DataFrame = pd.read_csv(TUNNEL_LOG)
        index  = df[ df['pdbid']==pdbid ].index.values.astype(int)
        if not index.size >0: 
            #doestn exist
            index = len(df.pdbid) + 1
            df.loc[index] = [pdbid, species, chosentunnels]
        else:
            # df.loc[index[0]] = [pdbid, species,chosentunnels]
            df.append([pdbid, species,chosentunnels]) 
        df.to_csv(TUNNEL_LOG,  index=False)
        csv_number = open(TUNNEL_RESULTS).read().split(',')
        print("Wrote {} -> {}, {} to {}".format(csv_number, pdbid,species, TUNNEL_LOG))
        tunneln.extend(csv_number)

    if tunneln[0]  in ['0', '-1']:
        print("Got 0 or -1, exiting.")
        exit(1)



def _lg_update(pdbid, **kwargs):
    """
    """





# The curated ids of tunnels(generated manually)
tunneln        = []






# Goal NOW: plot tunnel from the log.
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

