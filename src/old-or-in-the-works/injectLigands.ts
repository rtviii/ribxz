// import axios from "axios";
// import { Ligand, RibosomeStructure } from "../types";
// import { openMain, updateMain } from "../helpers/help";

// interface Response {
//   id: string;
//   ligandInfo: Array<PDBLigandinfo>;
// }

// interface PDBLigandinfo {
//   structureId : string;
//   chemicalID  : string;
//   chemicalName: string;
//   formula     : string;
//   smiles      : string;
// }
// export const injectLigands = async (pdbid: string) => {
//   const endpoint = `https://www.rcsb.org/pdb/json/ligandInfo?structureId=${pdbid}`;
//   const struct = openMain(pdbid);
//   var pdbresp: Response = await axios
//     .get(endpoint)
//     .then(response => response.data);

//   var ligandsreshaped: Array<Ligand> = pdbresp.ligandInfo.map(
//     (lig: PDBLigandinfo) => {
//       var newligand: Ligand = {
//         formula_weight: lig.formula,
//         pdbx_description: lig.
//         chemicalId   : lig.chemicalID,
//         chemicalName : lig.chemicalName,
//         cif_residueId: 'none',
//       };
//       return newligand;
//     }
//   );
//   const struct_ammended: RibosomeStructure = {
//     ...struct,
//     ...{ ligands: ligandsreshaped },
//   };
//   console.log("---- Inject Ligands: ");
//   updateMain(pdbid, struct_ammended)
// };
