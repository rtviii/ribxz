import axios from "axios";
import { chain, toString, uniq } from "lodash";
import {
  BanClass,
  Ligand,
  RibosomalProtein,
  RibosomeStructure,
  rRNA,
} from "../RibosomeTypes";
import { large_subunit_map } from "../resources/cumulative/large-subunit-map";
import { small_subunit_map } from "../resources/cumulative/small-subunit-map";
import {
  Nonpolymer_Entity,
  PDBGQLResponse,
  Polymer_Entity,
} from "./rcsb.gql.response-shape";

const getQuery = (pdbid: string) => {
  pdbid = pdbid.toUpperCase();
  return `https://data.rcsb.org/graphql?query={
    entry(entry_id: "${pdbid}") {
    rcsb_id
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_external_references {
      link
      type
      id
    }
    exptl {
      method
    }
    diffrn_source {
      details
      pdbx_synchrotron_beamline
      pdbx_synchrotron_site
      pdbx_wavelength
      pdbx_wavelength_list
      source
      type
    }
    em_3d_reconstruction {
      details
      algorithm
      resolution_method
      actual_pixel_size
      id
      resolution
      nominal_pixel_size
      refinement_type
      num_particles
      magnification_calibration
    }
    citation {
      rcsb_authors
      year
      title
      pdbx_database_id_DOI
    }
    struct_keywords {
      pdbx_keywords
      text
    }
    polymer_entities {
      entry {
        rcsb_id
      }
      pfams {
        rcsb_pfam_accession
        rcsb_pfam_comment
        rcsb_pfam_description
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        scientific_name
      }
      uniprots {
        rcsb_id
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
    }
    nonpolymer_entities {
      pdbx_entity_nonpoly {
        comp_id
        name
        entity_id
      }
      rcsb_nonpolymer_entity {
        formula_weight
        pdbx_description
      }
    }
  }
}`;
};

const matchPolymerNomenclature = (
  polymer: Polymer_Entity
): [BanClass[], string[]] | [[], []] => {
  if (!polymer.pfams) {
    return [[], []];
  } else {
    //   lookup in the fam-nom table
    var pfamIds = polymer.pfams.reduce<string[]>((acc, item) => {
      acc.push(item.rcsb_pfam_accession);
      return acc;
    }, []);

    var nomenclature: BanClass[] = [];
    pfamIds!.map(id => {
      Object.entries(large_subunit_map).map(entry => {
        if (entry[1].pfamDomainAccession.includes(id)) {
          nomenclature.push(entry[0] as BanClass);
        }
      });

      Object.entries(small_subunit_map).map(entry => {
        if (entry[1].pfamDomainAccession.includes(id)) {
          nomenclature.push(entry[0] as BanClass);
        }
      });
    });
    return [uniq(nomenclature), uniq(pfamIds)];
  }
};

const reshape_ToLigand = (nonpoly: Nonpolymer_Entity): Ligand => {
  return {
    pdbx_description: nonpoly.rcsb_nonpolymer_entity.pdbx_description,
    formula_weight: nonpoly.rcsb_nonpolymer_entity.formula_weight,
    chemicalId: nonpoly.pdbx_entity_nonpoly.comp_id,
    chemicalName: nonpoly.pdbx_entity_nonpoly.name,
    cif_residueId: "none",
  };
};

const reshape_ToRibosomalProtein = (
  polymers: Polymer_Entity[]
): RibosomalProtein[] => {
  return polymers.reduce(
    (rp_transformed: RibosomalProtein[], current_poly: Polymer_Entity) => {
      var [nomenclature, pfamIds] = matchPolymerNomenclature(current_poly);

      if (current_poly.pfams) {
        var pfam_comments = uniq(
          current_poly.pfams.map(pfam => pfam.rcsb_pfam_comment)
        );
        var pfam_descriptions = uniq(
          current_poly.pfams.map(pfam => pfam.rcsb_pfam_description)
        );
        var pfam_ids = uniq(
          current_poly.pfams.map(pfam => pfam.rcsb_pfam_accession)
        );
      } else {
        var pfam_comments: string[] = [];
        var pfam_descriptions: string[] = [];
        var pfam_ids: string[] = [];
      }

      var organism_ids: number[];
      var organism_descriptions: string[];
      organism_ids = current_poly.rcsb_entity_source_organism.map(
        org => org.ncbi_taxonomy_id
      );
      organism_descriptions = current_poly.rcsb_entity_source_organism.map(
        org => org.scientific_name
      );

      var uni = current_poly.uniprots
        ? current_poly.uniprots.map(entry => entry.rcsb_id)
        : [];

      // If the chain is duplicated, as it would be for two models packed inside a cif file from an xray experiment,
      // reduce each name to an individual record.
      var chains_in_record: RibosomalProtein[] = current_poly.entity_poly.pdbx_strand_id.includes(
        ","
      )
        ? current_poly.entity_poly.pdbx_strand_id.split(",").reduce(
            (chainobjects: RibosomalProtein[], chainid: string) => [
              ...chainobjects,
              {
                parent_rcsb_id: current_poly.entry.rcsb_id,
                pfam_comments: pfam_comments,
                pfam_descriptions: pfam_descriptions,
                pfam_accessions: pfam_ids,
                rcsb_source_organism_description: organism_descriptions,
                rcsb_source_organism_id: organism_ids,
                uniprot_accession: uni,
                rcsb_pdbx_description: current_poly.rcsb_polymer_entity
                  .pdbx_description
                  ? current_poly.rcsb_polymer_entity.pdbx_description
                  : null,
                entity_poly_strand_id: chainid,
                entity_poly_seq_one_letter_code:
                  current_poly.entity_poly.pdbx_seq_one_letter_code,
                entity_poly_seq_one_letter_code_can:
                  current_poly.entity_poly.pdbx_seq_one_letter_code_can,
                entity_poly_seq_length:
                  current_poly.entity_poly.rcsb_sample_sequence_length,
                entity_poly_polymer_type:
                  current_poly.entity_poly.rcsb_entity_polymer_type,
                entity_poly_entity_type: current_poly.entity_poly.type,
                nomenclature: nomenclature,
                surface_ratio: null,
              },
            ],
            []
          )
        : [
            {
              parent_rcsb_id: current_poly.entry.rcsb_id,
              pfam_comments: pfam_comments,
              pfam_descriptions: pfam_descriptions,
              pfam_accessions: pfam_ids,
              rcsb_source_organism_description: organism_descriptions,
              rcsb_source_organism_id: organism_ids,
              uniprot_accession: uni,
              rcsb_pdbx_description: current_poly.rcsb_polymer_entity
                .pdbx_description
                ? current_poly.rcsb_polymer_entity.pdbx_description
                : null,
              entity_poly_strand_id: current_poly.entity_poly.pdbx_strand_id,
              entity_poly_seq_one_letter_code:
                current_poly.entity_poly.pdbx_seq_one_letter_code,
              entity_poly_seq_one_letter_code_can:
                current_poly.entity_poly.pdbx_seq_one_letter_code_can,
              entity_poly_seq_length:
                current_poly.entity_poly.rcsb_sample_sequence_length,
              entity_poly_polymer_type:
                current_poly.entity_poly.rcsb_entity_polymer_type,
              entity_poly_entity_type: current_poly.entity_poly.type,
              nomenclature: nomenclature,
              surface_ratio: null,
            },
          ];
      return [...rp_transformed, ...chains_in_record];
    },
    []
  );
};
const reshape_TorRNA = (polymers: Polymer_Entity[]): rRNA[] => {
  return polymers.reduce(
    (rna_transformed: rRNA[], curr_poly: Polymer_Entity) => {
      var organism_ids: number[];
      var organism_descriptions: string[];

      if (curr_poly.rcsb_entity_source_organism) {
        organism_ids = curr_poly.rcsb_entity_source_organism.map(
          org => org.ncbi_taxonomy_id
        );
        organism_descriptions = curr_poly.rcsb_entity_source_organism.map(
          org => org.scientific_name
        );
      } else {
        organism_ids = [];
        organism_descriptions = [];
      }
      var chains_in_record: rRNA[] = curr_poly.entity_poly.pdbx_strand_id.includes(
        ","
      )
        ? curr_poly.entity_poly.pdbx_strand_id
            .split(",")
            .reduce((chains_in_record: rRNA[], chainid: string) => {
              return [
                ...chains_in_record,
                {
                  parent_rcsb_id: curr_poly.entry.rcsb_id,
                  rcsb_pdbx_description: curr_poly.rcsb_polymer_entity
                    .pdbx_description
                    ? curr_poly.rcsb_polymer_entity.pdbx_description
                    : null,
                  rcsb_source_organism_description: organism_descriptions,
                  rcsb_source_organism_id: organism_ids,
                  entity_poly_strand_id: chainid,
                  entity_poly_seq_one_letter_code:
                    curr_poly.entity_poly.pdbx_seq_one_letter_code,
                  entity_poly_seq_one_letter_code_can:
                    curr_poly.entity_poly.pdbx_seq_one_letter_code_can,
                  entity_poly_seq_length:
                    curr_poly.entity_poly.rcsb_sample_sequence_length,
                  entity_poly_polymer_type:
                    curr_poly.entity_poly.rcsb_entity_polymer_type,
                  entity_poly_entity_type: curr_poly.entity_poly.type,
                },
              ];
            }, [])
        : [
            {
              parent_rcsb_id: curr_poly.entry.rcsb_id,
              rcsb_pdbx_description: curr_poly.rcsb_polymer_entity
                .pdbx_description
                ? curr_poly.rcsb_polymer_entity.pdbx_description
                : null,
              rcsb_source_organism_description: organism_descriptions,
              rcsb_source_organism_id: organism_ids,
              entity_poly_strand_id: curr_poly.entity_poly.pdbx_strand_id,
              entity_poly_seq_one_letter_code:
                curr_poly.entity_poly.pdbx_seq_one_letter_code,
              entity_poly_seq_one_letter_code_can:
                curr_poly.entity_poly.pdbx_seq_one_letter_code_can,
              entity_poly_seq_length:
                curr_poly.entity_poly.rcsb_sample_sequence_length,
              entity_poly_polymer_type:
                curr_poly.entity_poly.rcsb_entity_polymer_type,
              entity_poly_entity_type: curr_poly.entity_poly.type,
            },
          ];

      return [...rna_transformed, ...chains_in_record];
    },
    []
  );
};

const inferOrganismFromPolymers: (
  proteins: Array<RibosomalProtein>
) => [string[], number[]] = proteins => {
  var structOrgIds: number[] = [];
  var structOrgNames: string[] = [];

  proteins.map(protein => {
    protein.rcsb_source_organism_description
      ? structOrgNames.push(...protein.rcsb_source_organism_description)
      : "";
    protein.rcsb_source_organism_id
      ? structOrgIds.push(...protein.rcsb_source_organism_id)
      : "";
  });

  return [uniq(structOrgNames), uniq(structOrgIds)];
};
const extractRefs = (
  external_refs: Array<{ link: string; type: string; id: string }>
) => {
  var externalRefIds: string[]   = [];
  var externalRefTypes: string[] = [];
  var externalRefLinks: string[] = [];

  external_refs.map(ref => {
    externalRefIds.push(ref.id);
    externalRefTypes.push(ref.type);
    externalRefLinks.push(ref.link);
  });
  return [externalRefIds, externalRefTypes, externalRefLinks];
};

export const requestGqlFrame = async (
  pdbid: string
): Promise<RibosomeStructure> => {
  var result = await axios.get(getQuery(pdbid)).then(response => {
    var pdbRecord: PDBGQLResponse = response.data.data;
    var proteins = pdbRecord.entry.polymer_entities.filter(
      poly => poly.entity_poly.rcsb_entity_polymer_type === "Protein"
    );
    var rnas = pdbRecord.entry.polymer_entities.filter(
      poly => poly.entity_poly.rcsb_entity_polymer_type === "RNA"
    );
    var ligands = pdbRecord.entry.nonpolymer_entities;

    var reshaped_proteins: RibosomalProtein[] = reshape_ToRibosomalProtein(
      proteins
    );
    var reshaped_rrnas: rRNA[] = reshape_TorRNA(rnas);
    var reshaped_ligands: Ligand[] | null =
      ligands == null ? null : ligands.map(r => reshape_ToLigand(r));

    var organismtuple = inferOrganismFromPolymers(reshaped_proteins);

    var organismNames = organismtuple[0];
    var organismIds = organismtuple[1];

    var externalRefs = extractRefs(pdbRecord.entry.rcsb_external_references);
    var em3d = pdbRecord.entry.em_3d_reconstruction
      ? pdbRecord.entry.em_3d_reconstruction[0]
      : null;
    var difn = pdbRecord.entry.diffrn_source
      ? pdbRecord.entry.diffrn_source[0]
      : null;

    var pub = pdbRecord.entry.citation[0];
    var kwords_text = pdbRecord.entry.struct_keywords
      ? pdbRecord.entry.struct_keywords.text
      : null;

    var reshaped: RibosomeStructure = {
      rcsb_id: pdbRecord.entry.rcsb_id,
      expMethod: pdbRecord.entry.exptl[0].method,
      resolution: pdbRecord.entry.rcsb_entry_info.resolution_combined[0],

      rcsb_external_ref_id: externalRefs[0],
      rcsb_external_ref_type: externalRefs[1],
      rcsb_external_ref_link: externalRefs[2],

      cryoem_exp_detail: em3d ? em3d.details : null,
      cryoem_exp_algorithm: em3d ? em3d.algorithm : null,
      cryoem_exp_resolution_method: em3d ? em3d.resolution_method : null,
      cryoem_exp_resolution: em3d ? em3d.resolution : null,
      cryoem_exp_num_particles: em3d ? em3d.num_particles : null,
      cryoem_exp_magnification_calibration: em3d
        ? em3d.magnification_calibration
        : null,

      diffrn_source_details: difn ? difn.details : null,
      diffrn_source_pdbx_synchrotron_beamline: difn
        ? difn.pdbx_synchrotron_beamline
        : null,
      diffrn_source_pdbx_synchrotron_site: difn
        ? difn.pdbx_synchrotron_site
        : null,
      diffrn_source_pdbx_wavelength: difn
        ? toString(difn.pdbx_wavelength)
        : null,
      diffrn_source_pdbx_wavelength_list: difn
        ? difn.pdbx_wavelength_list
        : null,
      diffrn_source_source: difn ? difn.source : null,
      diffrn_source_type: difn ? difn.type : null,

      citation_year: pub.year,
      citation_rcsb_authors: pub.rcsb_authors,
      citation_title: pub.title,
      citation_pdbx_doi: pub.pdbx_database_id_DOI,

      pdbx_keywords_text: kwords_text,

      _organismId: organismIds,
      _organismName: organismNames,

      proteins: reshaped_proteins,
      rnas: reshaped_rrnas,
      ligands: reshaped_ligands,
    };
    return reshaped;
  });
  return result;
};
