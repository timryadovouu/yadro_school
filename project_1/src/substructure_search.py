from rdkit import Chem
from fastapi import HTTPException

smiles_db = {
    "1": "c1ccccc1",
    "2": "CCO",
    "3": "CC(=O)Oc1ccccc1C(=O)O",
}

def substructure_search(substructure_smiles: str):
    substructure = Chem.MolFromSmiles(substructure_smiles)

    if substructure is None:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")
    
    matches = []
    for id in smiles_db:
        if Chem.MolFromSmiles(smiles_db[id]).HasSubstructMatch(substructure):
            matches.append({"id": id, "smiles": smiles_db[id]})
    return {"matched_molecules": matches}