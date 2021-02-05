import csv from "csv-parser";
import "shelljs";
import * as fs from "fs";
import path from "path";
import { RibosomalProtein, RibosomeStructure } from "../RibosomeTypes";

export const addSurfaceRatios = (pdbid: string) => {
  const results: any        = [];
  var   SURFACE_RATIOS_PATH = path.join(
    process.env.SURFACE_RATIOS as string,
    "surface_ratio_" + pdbid.toUpperCase() + ".csv"
  );
  if (!fs.existsSync(SURFACE_RATIOS_PATH)) {
    console.log(
      `addSurfaceRatios: ${SURFACE_RATIOS_PATH} does not exits. Skipping...`
    );
    return 0;
  }

  var currStruct: RibosomeStructure = JSON.parse(
    fs.readFileSync(
      path.join(
        process.env.STRUCTS_RCSB_GQL as string,
        pdbid.toUpperCase() + ".json"
      ),
      "utf8"
    )
  );
  fs.createReadStream(SURFACE_RATIOS_PATH)
    .pipe(csv())
    .on("data", (data: any) => results.push(data))
    .on("end", () => {
      //  For each row,  find the protein inside the record that matched the .name of the row
      // Assign the contents to the protein
      for (var row of results) {
        currStruct.proteins.map((protein: RibosomalProtein) => {
          if (protein.entity_poly_strand_id === row.name) {
            Object.assign(protein, { surface_ratio: parseFloat(row["ratio"]) });
          }
        });
      }
      // Write in-place
      fs.writeFileSync(
        path.join(
          process.env.STRUCTS_RCSB_GQL as string,
          pdbid.toUpperCase() + ".json"
        ),
        JSON.stringify(currStruct),
        "utf8"
      );
      console.log("Added surface ratios.");
    });
};
