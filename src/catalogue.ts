import fs from 'fs'
import path from 'path'
import { exit } from 'process'
import shell from 'shelljs'
import { addSurfaceRatios } from './rcsb-gql-api/addSurfaceRatios'
import { requestGqlFrame } from './rcsb-gql-api/requestGqlFrame'


const checkCatalogueExists = () =>{
    const catalogue = process.env.CATALOGUE as string
    if (!fs.existsSync(catalogue)){
        fs.closeSync(fs.openSync(catalogue, 'a'))
        fs.writeFileSync(catalogue, JSON.stringify({}))
        console.log("Struct catalogue is empty.");
        process.exit(0)
    }
}

interface CatalogueEntry{   
    [pdbid:string] :{
       reso   : number
       method : string
       ligands: string [] | null
       species: string[]
    }
}

export const restoreFromCatalogue = async () =>{
    const catalogue = process.env.CATALOGUE as string
    checkCatalogueExists()
    var allstructs = JSON.parse(fs.readFileSync(catalogue, 'utf-8'))
    for (var pdbid of Object.keys(allstructs)){
        console.log("Exists in the catalogue: " ,pdbid);
        await updateCatalogueWStruct(pdbid)
        await addSurfaceRatios(pdbid)
    }
}

export const seeCatalogue = () =>{
    const catalogue = process.env.CATALOGUE as string
    checkCatalogueExists()
    var allstructs = JSON.parse(fs.readFileSync(catalogue, 'utf-8'))
    console.log(allstructs);
}

export const updateCatalogueWStruct = async (pdbid: string) =>{

    const catalogue       = process.env.CATALOGUE as string
    const staticfiles     = process.env.STATIC_ROOT as string
    const target_filename = path.join(staticfiles,pdbid.toUpperCase(), pdbid.toUpperCase() + ".json")


    // Checking if the folder tree exists

    if (!fs.existsSync(path.dirname(target_filename))){
        console.log(path.dirname(target_filename), "does not exist. Creating folder tree.");
        shell.mkdir('-p', path.dirname(target_filename))
        
    }


    // if file already exists
    if (fs.existsSync(target_filename)){
        console.log(`Sturcture profile ${target_filename}.json EXISTS ALREADY.`);
        if (!fs.existsSync(catalogue)){
            fs.closeSync(fs.openSync(catalogue, 'a'))
            fs.writeFileSync(catalogue, JSON.stringify({}))
        }

    }
    // If not, grab it from rcsb
    else{

        var   rcsb_record_master = await requestGqlFrame(pdbid);
        fs.writeFileSync(target_filename, JSON.stringify(rcsb_record_master, null, 4))
        if (!fs.existsSync(catalogue)){
            fs.closeSync(fs.openSync(catalogue, 'a'))
            fs.writeFileSync(catalogue, JSON.stringify({}))
        }

        var allstructs = JSON.parse(fs.readFileSync(catalogue, 'utf-8'))

        var new_entry: CatalogueEntry = {
          [rcsb_record_master.rcsb_id.toUpperCase()]: {
            reso: rcsb_record_master.resolution,
            method: rcsb_record_master.expMethod,
            ligands: rcsb_record_master.ligands
              ? rcsb_record_master.ligands.map(( l:any ) => l.chemicalId)
              : [],
            species: rcsb_record_master._organismName,
          },
        };

        var allstructs = {...allstructs, ...new_entry}
        // Write to catalogue
        fs.writeFileSync(catalogue, JSON.stringify(allstructs, null, 4))
        console.log(`\n:=:-+-:===:@:===:+-+:=::=:-+-:===:@:===:+-+:=:=:-+-:===:@:===:+-+:=\nWrote ${rcsb_record_master.rcsb_id.toUpperCase()} to catalogue.`);
    }
}

