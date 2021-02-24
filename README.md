# ribxz 

Some tools and scripts to curate and analyze ribosomal structures as they pass from rcsb.org to ribosome.xyz.
TS cli in [ /src ](/src/cli.ts) mostly deals with pulling down and transforming the structures' representations from rcsb with own fields like the nomenclature. Some(most) information provided by  

Quoting from RCSB [ web services ](https://www.rcsb.org/docs/general-help/web-services-overview):
"All data stored in the PDB archive conform to the [PDBx/mmCIF data dictionary](https://mmcif.wwpdb.org/). This data is augmented with annotations coming from external resources and internally added fields. The RCSB PDB data representation, powered by the JSON Schema language, is connected to the [ data hierarchy ](https://data.rcsb.org/#data-organization)."

__The process of constructing an individual ribosome's profile is described in more detail [ here ](src/resources/StructureFromScratch.md).__



## Todos:

- Urgent, ordered:
    - [ ] Rewrite the schema in terms of pdbx/mmcif dictionary where possible, verify through 
    - [ ] Standardize tunnel extraction further around the [ logging utility ](/ciftools/TunnelScripts/TunnelLog.py). Some rewrites are due. 
    - [ ] Refine the definition of tunnel walls.
    - [ ] Better [ cataloguing ](catalogue_static_files.py) for static files.

- Later:

    - [ ] Centralize the pipeline for constructing a structure's profile along with the adjcent static files. From a PDBID to database, multiple transformations and script have to be applied, each lossy, which complicates the further down the line we go. Perhaps this can be a centralized log. Don't go crazy with this.
    - [ ] More robust mechanism for inducting a structure into the database. Currently, a poorly-parametrized [ bash script ](/src/resources/cypher-tools/induct_struct.sh).
    - [ ] Connect binding site, exit tunnel reports  to the db-schema such that they can be queried.


## Log:
