import axios from "axios";
// import { openMain, updateMain } from "../helpers/help";

interface Description {
    structureId  : string,
    title        : string,
    pubmedId     : string,
    expMethod    : string,
    organism     : string,
    replaces     : string,
    keywords     : string,
    nr_entities  : string,
    nr_residues  : string,
    nr_atoms     : string,
    publish_date : string,
    revision_date: string,
    status       : string,
}


export const injectMetadata = async (pdbid: string) => {
  const endpoint = `https://www.rcsb.org/pdb/json/describePDB?structureId=${pdbid}`;
  // const struct = openMain(pdbid);
  // var pdbresp:Array<Description> = await axios.get(endpoint).then(response => response.data);
  // var reshaped:RibosomeStructure = {...struct, ...{title: pdbresp[0].title, expMethod: pdbresp[0].expMethod} }
  
  // console.log("---- Inject Metadata: ");
  // updateMain(pdbid, reshaped)
};
