from fastapi import FastAPI, HTTPException, Depends
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError
from pydantic import BaseModel
from rdkit import Chem 
from models import Molecule
from database import engine, get_db
import os


app = FastAPI()

Molecule.metadata.create_all(bind=engine)

class MoleculeResponse(BaseModel):
    id: str
    smiles: str 

class MoleculeUpdate(BaseModel):
    smiles: str

@app.get("/")
def get_server():
    return {"server_id": os.getenv("SERVER_ID", "1")}

@app.post("/add")
def add_molecule(molecule: MoleculeResponse, db: Session = Depends(get_db)):
    if Chem.MolFromSmiles(molecule.smiles) is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES.")

    db_molecule = Molecule(id=molecule.id, smiles=molecule.smiles)
    db.add(db_molecule)
    try:
        db.commit()
    except IntegrityError:
        db.rollback()
        raise HTTPException(status_code=400, detail="Molecule already exists.")
    
    return {"message": f"Molecule with id {molecule.id} added successfully."}


@app.get("/search/{id}")
def get_molecule(id: str, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    
    return {"id": db_molecule.id, "smiles": db_molecule.smiles}


@app.get("/list")
def list_molecules(db: Session = Depends(get_db)):
    # return db.query(models.Molecule).all()
    return {"molecules": {molecule.id: molecule.smiles for molecule in db.query(Molecule).all()}}


@app.put("/update/{id}")
def update_molecule(id: str, molecule: MoleculeUpdate, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if db_molecule is None:
        raise HTTPException(status_code=400, detail="Invalid id")
    
    if Chem.MolFromSmiles(molecule.smiles) is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES.")
    
    if db_molecule.smiles == molecule.smiles:
        return {"message": f"Nothing changed for molecule with id: {id}"}

    db_molecule.smiles = molecule.smiles
    db.commit()
    return {"message": f"Successfully updated molecule with id: {id}. New value: {molecule.smiles}"}

@app.delete("/delete/{id}")
def delete_molecule(id: str, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if db_molecule is None:
        raise HTTPException(status_code=400, detail="Invalid id")

    db.delete(db_molecule)
    db.commit()
    return {"message": f"Successfully deleted molecule with id: {id}"}

@app.get("/substructure_search/{substructure_smiles}")
def substructure_search(substructure_smiles: str, db: Session = Depends(get_db)):
    substructure = Chem.MolFromSmiles(substructure_smiles)
    if substructure is None:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")
    
    matches = []
    for db_molecule in db.query(Molecule).all():
        if Chem.MolFromSmiles(db_molecule.smiles).HasSubstructMatch(substructure):
            matches.append({"id": db_molecule.id, "smiles": db_molecule.smiles})
    
    return {"matched_molecules": matches}