#!/usr/bin/env python3

import os,sys,numpy,json
from dotenv import load_dotenv

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    # Pass envfile-path to dotenv or other environ consumers.
    # envfile    = os.path.join(root(os.path.abspath(__file__),'ribxz'), '.env')
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')

import pandas as pd
from typing import List
from dotenv import load_dotenv
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from ciftools.structure import fetchStructure
from ciftools.neoget import _neoget

AAs          = ["ALA",'ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL']
Nucleotides  = ['A', 'U', 'G', 'C', 'T']
AMINO_ACIDS ={
"ALA"  :  0,
'ARG'  :  1,
'ASN'  :  0,
'ASP'  :  -1,
'CYS'  :  0,
'GLN'  :  0,
'GLU'  :  -1,
'GLY'  :  0,
'HIS'  :  0,
'ILE'  :  0,
'LEU'  :  0,
'LYS'  :  1,
'MET'  :  0,
'PHE'  :  0,
'PRO'  :  0,
'SER'  :  0,
'THR'  :  0,
'TRP'  :  0,
'TYR'  :  0,
'VAL'  :  0,
'SEC'  :  0,
'PYL'  :  0}

def get_CA_or(res:Residue)->Atom:

    """Returns the alpha carbon for the resiude if there is one, else the first atom in the list"""
    atoms:List[Atom] =list( res.get_atoms() )
    alphacarbons = list(filter(lambda atom: True if atom.get_name() == 'CA' else False ,atoms))
    return atoms[0] if len(alphacarbons) == 0 else alphacarbons[0]




class TunnelRecord:

    def __init__(self,
    data_path  : str,
    pdbid      : str,
    taxid      : int,
    csv_choices: List[int],
    comments   : List[str]=[""])->None: 

        """The record is initiated from csv files, but the walls are not rendered automatically."""

        self.struct_data_location   = data_path
        self.pdbid               = pdbid
        self.taxid               = taxid
        self.csv_choices         = csv_choices
        self.comments            = comments
        self.mole_dataframes:List[pd.DataFrame]     = []
        self.computed_properties = {"cumulaive_distance":  0}
        self.tunnel_walls:TunnelWalls        = {}
        self.has_rendered_walls = False

        distance          =  0

        mole_chosen_csvs  =  [os.path.join(self.struct_data_location,'csv', 'tunnel_{}.csv'.format(x)) for x in self.csv_choices]

        for f in mole_chosen_csvs:
            tunnel_instance   =  pd.read_csv(f)
            xyzr              =  tunnel_instance[['Distance','FreeRadius', 'X','Y','Z']]
            rows              =  len(xyzr)
            distance         +=  xyzr.iloc[rows-1]['Distance']
            self.mole_dataframes.append(tunnel_instance)
            
        self.computed_properties['cumulative_distance'] = distance

    def render_walls(self, radius:int):

        """Given a cif structure, render tunnel walls according to this record's  csv dataframes"""
        print("rendering walls")
        walls      = TunnelWalls(self.pdbid, fetchStructure(self.pdbid))

        total_data = pd.DataFrame()
        for f in self.mole_dataframes:
            total_data=total_data.append(f)

        walls.consumeMoleDataframe(total_data,radius)
        self.tunnel_walls = walls
        self.has_rendered_walls  = True

    def plot_radius(self)->None:
        df   = self.get_total_df()
        rd   = list(df['Radius'].array)
        dist = list(df['Distance'].array)
        plt.plot(dist, rd)

        plt.show()
    def get_total_df(self)->pd.DataFrame:
        """If multiple tunnels are present, concatenating into one, reversing radii."""

        if len(self.mole_dataframes) > 2:

            raise Exception("Too many tunnels found.")

        if len( self.mole_dataframes ) > 1:

            assert(len(self.mole_dataframes) == 2)

            pt1   = self.mole_dataframes[0][["Distance","Radius","X","Y","Z"]]
            length1 = pt1["Distance"][ len(pt1) -1] #Grabbing the total distance for first tunnel
            pt1['Distance'] = pt1[ 'Distance' ].values[::-1]
            pt1 = pt1[::-1]
            pt2   = self.mole_dataframes[1][["Distance","Radius","X","Y","Z"]]
            pt2['Distance']= pt2['Distance'].apply(lambda row: row +length1) # increase the length by first tunnel's; without reversing
            return pt1.append(pt2)

        else:
            assert(len(self.mole_dataframes)==1)
            return self.mole_dataframes[0][["Radius","X","Y","Z"]]

    def plot_all(self)->None:
        df          = self.get_total_df()
        xs          = df['X'].tolist()
        ys          = df['Y'].tolist()
        zs          = df['Z'].tolist()


        fig = plt.figure()
        ax  = plt.axes(projection='3d')
        ax.scatter3D(xs,ys,zs)
        ax.set_title(self.pdbid + ":" + ",".join(self.csv_choices))
        plt.show()


class TunnelWalls:

    def __init__(self, pdbid:str, structure: Structure) -> None:  
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
    def addResidue(self, res:Residue)->None:
        parentStrand = res.get_parent().get_id()
        if res.get_resname() not in [*AAs, *Nucleotides] and res not in self.ligands:
            self.ligands.append(res)

        if parentStrand not in self.adjacentRnaStrands and parentStrand not in self.adjacentRPStrands:

            response = _neoget("""
            match (n{{entity_poly_strand_id:"{parentStrand}"}})-[]-(r:RibosomeStructure{{rcsb_id:"{pdbid}"}}) \
            return {{type: n.entity_poly_polymer_type, nomenclature:n.nomenclature}};""".format_map({
                "parentStrand": parentStrand,
                "pdbid"      : self.pdbid
            }))
            try:
                profile = response[0]
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
                print("Got nbr ", nbr)
                self.addResidue(nbr)

        df.apply(getVicinity, axis=1)
        print("Total unique residues: ",self.rescount)

    def generateReport(self, write_to_path:str=""):
        """Consume Mole Dataframe first. Things are empty otherwise. Should be in the appropriate record."""

        def getResInfo(res:Residue, nomMap, polytype:str ): 

            parent       : Chain.Chain = res.get_parent()
            parentStrand: str          = parent.get_id()
            resid        : int         = res.get_id()[1]
            resname      : str         = res.get_resname()
            rescoord                   = get_CA_or(res).get_coord().tolist()
            

            nom   = nomMap[parentStrand] if polytype == "Protein" else None 
            islig = False if resname.upper() in [*Nucleotides, *AAs] else True
            print("Added residue ", resname)

            return {
                "strand"  : parentStrand,
                "resid"   : resid,
                "resname" : resname,
                "polytype": polytype,
                "nom"     : nom,
                "isligand": islig,
                "rescoord": rescoord
                }
        
        protwall  =  {}
        rnawall   =  {}
        presentLigands = [getResInfo(l,nomMap=self.nomenclatureMap, polytype="Other") for l in self.ligands]

        for tpl in self.getProteinResidues().items():
            protwall[ tpl[0] ] = [getResInfo(x,nomMap=self.nomenclatureMap, polytype="Protein") for x in tpl[1]]

        for tpl in self.getRnaResidues().items():
            rnawall[ tpl[0] ] = [getResInfo(x,nomMap=self.nomenclatureMap, polytype="RNA") for x in tpl[1]]

        report = {
            "pdbid"      : self.pdbid,
            "probeRadius": self.radius,
            "rna"        : rnawall,
            "proteins"   : protwall,
            "ligands"    : presentLigands,
            "nomMap"     : self.nomenclatureMap}

        if write_to_path != "":

            FILENAME  =  "{}_TUNNEL_REPORT.json".format(self.pdbid)
            OUTPATH   =  os.path.join(write_to_path,FILENAME)

            with open(OUTPATH, "w") as outfile:
                json.dump(report, outfile)
                print("Has written to path {}".format(OUTPATH))
        return report

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

        taxid = row['TaxId'].values[0]
        molechoices = row['MoleStatus'].values[0].split(',')

        return TunnelRecord(
            os.path.join(self.__tunnels_path,str(taxid),pdbid),
            pdbid,
            taxid,
            molechoices)


# PDBID=sys.argv[1].upper()
# PDBID='3JCD'
# LOG = Log(os.getenv('TUNNEL_LOG'))
# """You have to render walls on a given records first before you generated the report or do anything else."""
# record_3jcd = LOG._struct(PDBID)
# # x.plot_all()
# # 
# record_3jcd.render_walls(10)
# # STRUCT_REPORT_PATH = os.path.join(os.getenv('STATIC_ROOT'),PDBID,'TUNNEL', f'{PDBID}_TUNNEL_REPORT.json')
# x = record_3jcd.tunnel_walls.generateReport(write_to_path='.')



