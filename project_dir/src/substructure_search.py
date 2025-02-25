from rdkit import Chem
from fastapi import HTTPException
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

smiles_db = {
    "1": "c1ccccc1",
    "2": "CCO",
    "3": "CC(=O)Oc1ccccc1C(=O)O",
}

def substructure_search(substructure_smiles: str):
    logger.info("call substructure_search function")
    substructure = Chem.MolFromSmiles(substructure_smiles)

    if substructure is None:
        logger.error("invalid substructure SMILES")
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")
    
    matches = []
    for id in smiles_db:
        
        if Chem.MolFromSmiles(smiles_db[id]).HasSubstructMatch(substructure):
            logger.info(f"adding {smiles_db[id]} as match")
            matches.append({"id": id, "smiles": smiles_db[id]})
    return {"matched_molecules": matches}