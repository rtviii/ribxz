export interface Polymer_Entity {
  entry:{
    rcsb_id:string
  }
  pfams: Array<{
    rcsb_pfam_accession  : string;
    rcsb_pfam_comment    : string;
    rcsb_pfam_description: string;
    }> | null;

  rcsb_entity_source_organism: Array<{ 
    ncbi_taxonomy_id: number;
    scientific_name : string;
    }>;

  uniprots: Array<{rcsb_id: string;}> | null;

  rcsb_polymer_entity: {
    pdbx_description: string;
  };

  entity_poly: {
    pdbx_seq_one_letter_code    : string;
    pdbx_seq_one_letter_code_can: string;
    pdbx_strand_id              : string;
    rcsb_entity_polymer_type    : string;
    rcsb_sample_sequence_length : number;
    type                        : string;
  };
}
export interface Nonpolymer_Entity {
  pdbx_entity_nonpoly: {
    comp_id  : string;
    name     : string;
    entity_id: string;
  };
  rcsb_nonpolymer_entity: {
    formula_weight  : number;
    pdbx_description: string;
  };
}

export interface PDBGQLResponse {
  entry: {
    rcsb_id: string;
    rcsb_entry_info         : {resolution_combined:Array<number>;};
    rcsb_external_references: Array<{link:string; type:string; id:string}>
    exptl                   : Array<{ method: string }>;

    diffrn_source           : Array<{
      details                  : string | null;
      pdbx_synchrotron_beamline: string |null ;
      pdbx_synchrotron_site    : string | null;
      pdbx_wavelength          : string | number | null;
      pdbx_wavelength_list     : string | null;
      source                   : string | null;
      type                     : string | null;
    }>;

    em_3d_reconstruction: Array<{
      details                  : string | null;
      algorithm                : string | null;
      resolution_method        : string | null;
      actual_pixel_size        : number | null;
      id                       : string | number;
      resolution               : number | null;
      refinement_type          : null | string;
      num_particles            : number | null;
      magnification_calibration: null | string;
    }>;

    citation: Array<{
      rcsb_authors        : Array<string>;
      year                : number;
      title               : string;
      pdbx_database_id_DOI: string;
    }>;
    struct_keywords: {
      text         : string;
    };
    polymer_entities   : Array<Polymer_Entity>;
    nonpolymer_entities: Array<Nonpolymer_Entity> | null;
  };
}
