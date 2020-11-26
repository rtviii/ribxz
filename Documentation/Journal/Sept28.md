- Adding a function for pulling and integrating ligands from PDB into the schema
- Not sure currently what's the best way to capture the position of the ligand in the schema. Opting for the residueId as per Biopython in this case:

```typescript
export interface Ligand {
  name: string;
  chemicalId: string;
  neighbors: Array<string>;
  cif_residueId: number;
}
```

---
````javascript
// I want to get rid of commander already. Will rewrite the parser at some point to have custom-types for commands, arg-shapes for each...
````
---
Cleaning the malformed/missing/erroneous nomenclature up, uploading to neo4j.
(To detect where things went wrong : 

````jq
jq '.proteins[] | select(.nomenclature| length > 1)| {c:._PDBChainId,desc:.description,n:.nomenclature, name:._PDBName,}' 4UG0.json
```` 
---

__eL8__ nomenclature class discovered by inducion through 5T2A: Zhang et al. annotation follow Ban's nomenclature.