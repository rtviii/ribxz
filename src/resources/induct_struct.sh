#!/usr/bin/bash



NEOIMPORT='/var/lib/neo4j/import'
USER='rt'
PASS='55288'
# The current instance of the neo4j database.
DATABASE_NAME='ribban03'

#The assumption is that the static folder with all the json and cif files 
#organized by structure exists exactly at $NEOIMPORT/static/

#⋱⋰⋯⋯⋯⋅⋄⋱⋰⋯⋯⋯⋱⋰⋯⋯⋯⋅⋄⋱⋰⋯⋯⋯⋱⋰⋯⋯⋯⋅⋄⋱⋰⋯⋯⋯⋱⋰⋯⋯⋯⋅⋄⋱⋰⋯⋯⋯

filepath=$1

if [ -f $filepath ];
then
	file=$(basename $filepath)
	extension=${file: -4}
	if [ $extension != "json" ];
	then
		echo "The profile file must be a json registry. Exiting."
		exit 2
	fi
	# verify that the id is uppercase
	structid=${file::4}
	structid=${structid^^}

else
	echo "$filepath is not an acceptable file"
	exit -1
fi

#Throwing struct into cypher-shell
echo "call apoc.load.json(\"file:///static/$structid/$file\") yield value
	with    value.rcsb_id                                      as pdbid,
        value.expMethod                                    as exp, 
        value.resolution                                   as reso, 
        value.rcsb_external_ref_id                              as ref_id,
        value.rcsb_external_ref_type                            as ref_type,
        value.rcsb_external_ref_link                            as ref_link,
        value.cryoem_exp_detail                            as c_exp_detail,
        value.cryoem_exp_algorithm                          as c_exp_algo,
        value.cryoem_exp_resolution                         as c_exp_reso,
        value.cryoem_exp_resolution_method                  as c_exp_reso_method,
        value.cryoem_exp_magnification_calibration          as c_exp_mag_calibration,
        value.cryoem_exp_num_particles                      as c_exp_n_particles,
        value.diffrn_source_details                   as difn_details,
        value.diffrn_source_pdbx_synchrotron_beamline as difn_sync_beam,
        value.diffrn_source_pdbx_synchrotron_site     as difn_sync_site,
        value.diffrn_source_pdbx_wavelength           as difn_omega,
        value.diffrn_source_pdbx_wavelength_list      as difn_omega_list,
        value.diffrn_source_source                    as difn_source,
        value.diffrn_source_type                      as difn_type,
        value.citation_year                            as cit_year,
        value.citation_rcsb_authors                    as cit_authors,
        value.citation_title                           as cit_title,
        value.citation_pdbx_doi                        as cit_doi,
        value.pdbx_keywords_text as kwords, 
        value._organismId as orgid, 
        value._organismName as orgname, 
        value
merge (struct:RibosomeStructure{
        rcsb_id                                : pdbid,
        expMethod                              : exp,
        resolution                             : reso,

        citation_year                          : cit_year,
        citation_rcsb_authors                  : cit_authors,
        citation_title                         : cit_title,
        citation_pdbx_doi                      : cit_doi,

        _organismId                            : orgid,
        _organismName                          : orgname
        })
        on create set
        struct.cryoem_exp_detail                      = CASE WHEN c_exp_detail  = null then "null" else c_exp_detail END,
        struct.cryoem_exp_algorithm                   = CASE WHEN c_exp_algo  = null then "null" else c_exp_algo END,
        struct.cryoem_exp_resolution                  = CASE WHEN c_exp_reso  = null then "null" else c_exp_reso END,
        struct.cryoem_exp_resolution_method           = CASE WHEN c_exp_reso_method  = null then "null" else c_exp_reso_method END,
        struct.rcsb_external_ref_id                   = CASE WHEN ref_id = null then "null" else ref_id END,
        struct.rcsb_external_ref_type                 = CASE WHEN ref_type = null then "null" else ref_type END,
        struct.rcsb_external_ref_link                 = CASE WHEN ref_link = null then "null" else ref_link END,
        struct.cryoem_exp_magnification_calibration   = CASE WHEN c_exp_mag_calibration  = null then "null" else c_exp_mag_calibration END,
        struct.cryoem_exp_num_particles               = CASE WHEN c_exp_n_particles  = null then "null" else c_exp_n_particles END,
        struct.diffrn_source_details                  = CASE WHEN difn_details  = null then "null" else difn_details END,
        struct.diffrn_source_pdbx_synchrotron_beamline= CASE WHEN difn_sync_beam  = null then "null" else difn_sync_beam END,
        struct.diffrn_source_pdbx_synchrotron_site    = CASE WHEN difn_sync_site  = null then "null" else difn_sync_site END,
        struct.diffrn_source_pdbx_wavelength          = CASE WHEN difn_omega  = null then "null" else difn_omega END,
        struct.diffrn_source_pdbx_wavelength_list     = CASE WHEN difn_omega_list  = null then "null" else difn_omega_list END,
        struct.diffrn_source_source                   = CASE WHEN difn_source  = null then "null" else difn_source END,
        struct.diffrn_source_type                     = CASE WHEN difn_type  = null then "null" else difn_type END
with value, struct
unwind value.proteins as protein
with    protein,
        value,
        struct
merge (rp:RibosomalProtein {
        parent_rcsb_id                     : protein.parent_rcsb_id,
        pfam_comments                      : protein.pfam_comments,
        pfam_descriptions                  : protein.pfam_descriptions,
        pfam_accessions                    : protein.pfam_accessions,
        rcsb_source_organism_description   : protein.rcsb_source_organism_description,
        rcsb_source_organism_id            : protein.rcsb_source_organism_id,
        uniprot_accession                  : protein.uniprot_accession,
        rcsb_pdbx_description              : protein.rcsb_pdbx_description,
        entity_poly_strand_id              : protein.entity_poly_strand_id,
        entity_poly_seq_one_letter_code    : protein.entity_poly_seq_one_letter_code,
        entity_poly_seq_one_letter_code_can: protein.entity_poly_seq_one_letter_code_can,
        entity_poly_seq_length             : protein.entity_poly_seq_length,
        entity_poly_polymer_type           : protein.entity_poly_polymer_type,
        entity_poly_entity_type            : protein.entity_poly_entity_type,
        nomenclature                       : protein.nomenclature
    })-[:RibosomalProtein_of]->(struct)
    on create set 
    rp.rcsb_pdbx_description = CASE WHEN protein.rcsb_pdbx_description = null then "null" else protein.rcsb_pdbx_description END
with rp, struct, value
unwind rp.pfam_accessions as pfamils
match (pf:PFAMFamily {family_id:pfamils})
with rp,struct,value,pf
merge (rp)-[:Belogns_To]->(pf)
with value,struct
unwind value.rnas as rna
merge (s:rRNA{
        parent_rcsb_id                     :rna.parent_rcsb_id,
        rcsb_source_organism_description   :rna.rcsb_source_organism_description,
        rcsb_source_organism_id            :rna.rcsb_source_organism_id,
        entity_poly_strand_id              :rna.entity_poly_strand_id,
        entity_poly_seq_one_letter_code    :rna.entity_poly_seq_one_letter_code,
        entity_poly_seq_one_letter_code_can:rna.entity_poly_seq_one_letter_code_can,
        entity_poly_seq_length             :rna.entity_poly_seq_length,
        entity_poly_polymer_type           :rna.entity_poly_polymer_type,
        entity_poly_entity_type            :rna.entity_poly_entity_type
        })-[:rRNA_of]->(struct)
on create set s.rcsb_pdbx_description = CASE WHEN rna.rcsb_pdbx_description = null then "null"  else rna.rcsb_pdbx_description END
with value, struct
unwind value.ligands as lig
merge (l:Ligand {
        chemicalId      : lig.chemicalId,
        chemicalName    : lig.chemicalName,
        cif_residueId   : lig.cif_residueId,
        formula_weight  : lig.formula_weight,
        pdbx_description: lig.pdbx_description})
merge (l)<-[:ContainsLigand{cif_residueId: lig.cif_residueId}]-(struct)" | cypher-shell -u $USER -p $PASS -d $DATABASE_NAME --format plain
