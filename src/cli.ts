import program from "commander";
// import { inspect } from "./";
import { restoreFromCatalogue, seeCatalogue, updateCatalogueWStruct } from "./catalogue";

const populateCLI = () => {
  program.description(" ▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯  RIBXYZ ON  ▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯ \n");
  // this is named assistance because the help function is already defined on the program api.
  program.option("-h, --assistance");
  program.option("-gist, --displayGist");
  program.option("-struct, --generateStructureProfile [pdbid]");
  program.option("-cat, --catalogue");

// Not implemented
  program.option("-catrestore, --restoreFromCatalogue");
  program.option("-catdownload, --downloadStaticCifsCatalogue");
};

export const cli = (cargs: Array<string>) => {
  populateCLI();
  program.parse(cargs);
  ribxz_exec();
};

const gist = ()=>{

  console.log("ORGANIZATION\n" + 
  "The main storage directory for structs is always [static] at the root of the application.\n"+
  "Each structure's folder is strictly in UPPERCASE and is its 4-letter RCSB code.\n" + 
  "Within the struct folder XXXX.json structure profile is to be found along with the XXXX.cif crystallographic file.\n" +
  "The structur folder might also include CHAINS and LIGANDS subfolders that contain these data parsed with Biopython.\n")

  console.log("\n------------------ENVIRONMENT---------------")
  console.log("Should be defined in the .env file in the root dir. Relevant variables:")
  console.log(`[ STATIC_ROOT ]: now set to ${process.env.STATIC_ROOT}\n The main repository for static files and struct profiles.`)
  console.log(`[ CATALOGUE ]: now set to ${process.env.CATALOGUE}\n The catalogue file containing structs currently profiled. `)
  console.log(`[ CYPHER_IMAGE ]: now set to ${process.env.CYPHER_IMAGE}\n The Cypher script corresponding to the "types.ts" definition of a structure and sufficient to commit the profil to Neo4j in one go. `)
  console.log("--------------------------------------------\n")
  
  console.log("### Mechanism ###\n")
  console.log("The moving pieces are thus:")
  console.log("I pull the structure data from RCSB GraphQL API and reshape it according to the schema defined in `types.ts`. This is the authoriative file on how the object should be formed.")
  console.log("As a part of profile-building, i commit a short summary of the structure to the catalogue." )
  console.log("There is a phase of splitting the structs into chains and ligands which relies on the python scripts but should end up in the STATIC eventually.")
  console.log("The static should be then moved to the appropirate path on the server to serve the .cif")
  console.log("The profiles are inducted to the database by the `induct_struct.sh` script in Resources.")

  console.log("\n\nIf you are not finding something here, take a look at the Notes.md...")
  console.log("A convenient way to get all the structs is by generating a custom report on RCSB's search and then\n AWK'ing the struct names.")
}

const ribxz_exec = async () => {


  if (program.assistance){
    console.table({
      help: { invocation: "-h, --assistance", effect: "", comments: "" },
      struct: {
        invocation: "-struct, --generateStructureProfile [pdbid]",
        effect:
          'Pull a structure from RCSB PDB according to the "scheme" defined in the /types.ts',
        comments: "Relevant variables: ",
      },
      catalogue: {
        invocation: "-cata, --catalgoue",
        effect: "Display short-format catalogue",
        comments: "1. Should probably have some sort of anchor in a csv file of just struct names\n 2.Restore from that ",
      },
    });
  }

  if (program.generateStructureProfile) {
    var pdbid:string = program.generateStructureProfile
    if (pdbid.length < 4){
      console.log("Enter a valid RCSB PDB ID to build from RCSB's gql api.");
      process.exit(2)
    }
    await updateCatalogueWStruct(pdbid)
  }

  if (program.catalogue){
    seeCatalogue()
  }
};
