import pytest
from src.substructure_search import substructure_search
from fastapi import HTTPException


def test_valid_substructure_search():
    substructure_smiles = "c1ccccc1" 
    result = substructure_search(substructure_smiles)
    assert len(result["matched_molecules"]) == 2
    assert result["matched_molecules"][0]["id"] == "1"

def test_no_match_substructure():
    substructure_smiles = "CCF"  
    result = substructure_search(substructure_smiles)
    assert len(result["matched_molecules"]) == 0

def test_multiple_matches():
    substructure_smiles = "CC"  
    result = substructure_search(substructure_smiles)
    assert len(result["matched_molecules"]) == 2
    assert set(mol["id"] for mol in result["matched_molecules"]) == {"2", "3"}

def test_invalid_smiles_input():
    substructure_smiles = "invalid_smiles"
    with pytest.raises(HTTPException) as excinfo:
        substructure_search(substructure_smiles)
    assert excinfo.value.status_code == 400
    assert excinfo.value.detail == "Invalid substructure SMILES."