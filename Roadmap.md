


## Pre-Processing of Structures  (Urgent)

[ __Guide__ ](https://github.com/rtviii/ribxz/blob/master/src/resources/StructureFromScratch.md)

- [ ] A single script for injesting new structures into the database. 
    
    Though mostly automated, it still currently takes a few steps to pick the structures from PDB, transform them and generate the relevant data. This needs to be pushed to a point where a cronjob that runs twice a month is sufficient.
    (**CI**). This is problematic at the following stages:

    - I have been frequently altering/adding fields to the schema as new traits were added. This should not be a problem after the .cif dictionary is adopted.

    - Missing fields and irregular structures come up fairly frequently and the preprocessing fails silently then.  Have to be addressed not to accumulate errors downstream.

    - Not all components nominally included in the structure's profile are extractable(mainly an issue with the parser, again, and certain ligands' binding sites). Hence, static files have to be logged to have a handle on what can actually be downloaded.


# Frontend

### General

- [ ] [ ? depends on the parser ] Viewer Integration  ( __High priority__ )
- [ ] [ 1-2 weeks ]               Deep search, filtering logic  ( __High priority__ )
- [ ] [ 1 week ]                  Corrections to layout   ( __High priority__ )
- [ ] [ 1 day ]                   Tutorial               ( __High priority__ )

### Protein Profiling into a Workspace Page 

- [ ] [ 3-5 days ] A unified interface (Proteovision-like) on a single structure with metrics and data on each subcomponent.

### Detected Issues:

+ Ligand/Structure page, nbrs representation doesn't always switch to the next struct. Probably selection errors in the front


# Backend


- [ ] [ 2-3 days ] A method for constructing a valid pdbx/mmcif structure either via [ Bio.PDB.Fast ](https://biopython.org/docs/1.75/api/Bio.PDB.MMCIFParser.html) or [ gemmi ](https://gemmi.readthedocs.io/en/latest/). ( __High priority__ )

    Currently the parsing/packing is done via PyMol's cli and the structures saved by it are ill-formed (mising the regular .mmcif structure header). This is a problem that has been of interest for a while now both for us and Anton's group and what the subcomponent alignment feature depends on.

- [ ] [ 2 day ] Rewrite ribosome types in terms of the .cif-dictionary terms where possible. Relatively small adjustment that would go a long way towards broader usage of RCSB's gql api at the front and back, eliminate schema disputes at the pre-processing step. *When in doubt, stick with their [ data model ](https://data.rcsb.org/index.html#data-organization)* 

- [ ] [ 0.5 day ] Bulk downloads (zip, csv). (__High priority__) 
- [ ] [ 0.5 day ] Additional formats (csv, fasta, fastaq,... )

- Integrations with Proteovision (__Incremental__) 


# Data Extraction


- [ ] [ 0.5 day] Integrating Tunnel Log into the database. (__Incremental__)
- [ ] [ 1-2 days ]Integrating Ligand Neighborhoods into the database  (__Incremental__)


----------- 

For Cole: 

- Translating the Matlab script (__Incremental__)
- Ligand-parsing script (__Incremental__)
- Occasional Reviews/sweeps for bugs, second perspective
- Bioinformatics stack exchange, favorable formats 