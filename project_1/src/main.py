from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem 
from os import getenv

app = FastAPI()
smiles_db = {
    "1": "c1ccccc1",
    "2": "CCO",
    "3": "CC(=O)Oc1ccccc1C(=O)O",
  }

class Molecule(BaseModel):
    id: str
    smiles: str

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

@app.post("/add")
def add_molecule(molecule: Molecule):
    if molecule.id in smiles_db:
        raise HTTPException(status_code=400, detail="Molecule already exists.")
    
    if Chem.MolFromSmiles(molecule.smiles) is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES.")

    smiles_db[molecule.id] = molecule.smiles
    return {"message": f"molecule with id {molecule.id} added successfully."}

@app.get("/search/{id}")
def get_molecule(id: str):
    smiles = smiles_db.get(id)
    if smiles is None:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    return {"id": id, "smiles": smiles}

@app.get("/list")
def list_molecules():
    return {"molecules": smiles_db}

@app.put("/update/{id}")
def update_molecule(id: str, molecule: Molecule):
    if id not in smiles_db:
        raise HTTPException(status_code=400, detail="Invalid id")
    
    if Chem.MolFromSmiles(molecule.smiles) is None:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")
    
    smiles_db[id] = molecule.smiles
    return {"message": f"sucess update for molecule with id: {id}. new value: {molecule.smiles}"}

@app.delete("/delete/{id}")
def delete_molecule(id: str):
    if id not in smiles_db or None:
        raise HTTPException(status_code=400, detail="Invalid id")
    
    del smiles_db[id]
    return {"message": f"sucess delete molecule with id: {id}"}

@app.get("/substructure_search/{substructure_smiles}")
def substructure_search(substructure_smiles: str):
    substructure = Chem.MolFromSmiles(substructure_smiles)

    if substructure is None:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")
    
    matches = []
    for id in smiles_db:
        if Chem.MolFromSmiles(smiles_db[id]).HasSubstructMatch(substructure):
            matches.append({"id": id, "smiles": smiles_db[id]})
    return {"matched_molecules": matches}