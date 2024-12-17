from celery_worker import celery 
from rdkit import Chem
from sqlalchemy.orm import Session
from models import Molecule, TaskResult
from database import get_db, SessionLocal
from fastapi import Depends, HTTPException
import json
import time

@celery.task
def substructure_search_task(substructure_smiles: str):
    substructure = Chem.MolFromSmiles(substructure_smiles)
    if substructure is None:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")
    
    cache_key = f"sub_search:{substructure_smiles}" 
    cached_result = celery.backend.get(cache_key)
    if cached_result:
        result = json.loads(cached_result) 
        return {"source": "cache", **result}
    
    db = SessionLocal()
    try:
        time.sleep(10)
        matches = []
        for db_molecule in db.query(Molecule).all():
            if Chem.MolFromSmiles(db_molecule.smiles).HasSubstructMatch(substructure):
                matches.append({"id": db_molecule.id, "smiles": db_molecule.smiles})
        response_data = {"matched_molecules": matches}
        celery.backend.set(cache_key, json.dumps(response_data)) 
        return {"source": "database", "matched_molecules": matches}
    finally:
        db.close()
