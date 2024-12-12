from fastapi import FastAPI, HTTPException, Depends
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from pydantic import BaseModel
from rdkit import Chem 
from models import Molecule
from database import engine, get_db
import os
import logging
import redis 
import json

app = FastAPI()

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

Molecule.metadata.create_all(bind=engine)

# connect to Redis
redis_client = redis.Redis(host='redis', port=6379, db=0)

def get_cached_result(key: str): 
    result = redis_client.get(key) 
    if result: 
        return json.loads(result) 
    return None 
 
# refresh time in seconds (default: 1 minute)
def set_cache(key: str, value: dict, expiration: int = 60): 
    redis_client.setex(key, expiration, json.dumps(value))

class MoleculeResponse(BaseModel):
    id: str
    smiles: str 

class MoleculeUpdate(BaseModel):
    smiles: str

@app.get("/")
def get_server():
    logger.info("getting server id")
    return {"server_id": os.getenv("SERVER_ID", "1")}

@app.post("/add")
def add_molecule(molecule: MoleculeResponse, db: Session = Depends(get_db)):
    if Chem.MolFromSmiles(molecule.smiles) is None:
        logger.error(f"invalid SMILES: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES.")

    db_molecule = Molecule(id=molecule.id, smiles=molecule.smiles)
    db.add(db_molecule)
    logger.info(f"trying to add molecule {molecule.smiles} with id: {molecule.id} to database")
    try:
        db.commit()
    except IntegrityError:
        db.rollback()
        logger.error(f"molecule {molecule.smiles} already exists")
        raise HTTPException(status_code=400, detail="Molecule already exists.")
    
    logger.info(f"successfull added molecule {molecule.smiles} with id: {molecule.id}")
    return {"message": f"Molecule with id {molecule.id} added successfully."}


@app.get("/search/{id}")
def get_molecule(id: str, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if db_molecule is None:
        logger.error("molecule not found")
        raise HTTPException(status_code=404, detail="Molecule not found.")
    logger.info(f"successfull search for molecule with id: {db_molecule.id} and smiles: {db_molecule.smiles}")
    return {"id": db_molecule.id, "smiles": db_molecule.smiles}


@app.get("/list")
def list_molecules(limit: int = 100, offset: int = 0, db: Session = Depends(get_db)):
    logger.info(f"listing molecules with limit {limit} starting from offset {offset}")

    cache_key = f"list_molecules:{limit}:{offset}" 
    cached_result = get_cached_result(cache_key)
    if cached_result: 
        logger.info(f"getting data from cache: {cache_key}")
        return {"source": "cache", "data": cached_result}
    
    logger.info(f"getting data from database: {cache_key}")
    molecules_number = db.query(Molecule).count()

    def molecule_iterator(limit, offset):
        query = db.query(Molecule).offset(offset).limit(limit)
        for molecule in query: 
            # yield {"id": molecule.id, "smiles": molecule.smiles}
            yield molecule

    molecules = {molecule.id: molecule.smiles for molecule in molecule_iterator(limit, offset)}
    response_data = {"number_of_molecules": molecules_number, "molecules": molecules}
    set_cache(cache_key, response_data) 

    # return {"molecules": {molecule.id: molecule.smiles for molecule in db.query(Molecule).all()}}
    return {"source": "database", "data": response_data}

@app.put("/update/{id}")
def update_molecule(id: str, molecule: MoleculeUpdate, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if db_molecule is None:
        logger.error(f"invalid id")
        raise HTTPException(status_code=400, detail="Invalid id")
    
    if Chem.MolFromSmiles(molecule.smiles) is None:
        logger.error(f"invalid SMILES: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES.")
    
    if db_molecule.smiles == molecule.smiles:
        logger.error(f"was given the same molecule with id {id} (value: {molecule.smiles})")
        return {"message": f"Nothing changed for molecule with id: {id}"}

    db_molecule.smiles = molecule.smiles
    db.commit()
    logger.info(f"successfull update for molecule with id: {id}, new value: {molecule.smiles}")
    return {"message": f"Successfully updated molecule with id: {id}. New value: {molecule.smiles}"}

@app.delete("/delete/{id}")
def delete_molecule(id: str, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if db_molecule is None:
        logger.error(f"cannot delete molecule with id {id}")
        raise HTTPException(status_code=400, detail="Invalid id")

    db.delete(db_molecule)
    db.commit()
    logger.info(f"successfull delete molecule with id {id}")
    return {"message": f"Successfully deleted molecule with id: {id}"}

@app.get("/substructure_search/{substructure_smiles}")
def substructure_search(substructure_smiles: str, db: Session = Depends(get_db)):
    substructure = Chem.MolFromSmiles(substructure_smiles)
    if substructure is None:
        logger.error("invalid substructure SMILES")
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")
    
    matches = []
    for db_molecule in db.query(Molecule).all():
        if Chem.MolFromSmiles(db_molecule.smiles).HasSubstructMatch(substructure):
            logger.info(f"adding {db_molecule.smiles} as match for {substructure_smiles}")
            matches.append({"id": db_molecule.id, "smiles": db_molecule.smiles})
    
    return {"matched_molecules": matches}
