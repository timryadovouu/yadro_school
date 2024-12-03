from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError
import pytest
import os

from src.main import app
from src import models

TEST_DATABASE_URL = "postgresql://${POSTGRES_USER}:${POSTGRES_PASSWORD}@db:5432/${POSTGRES_DB}"
# TEST_DATABASE_URL = "postgresql://fastapi_newbie:password@localhost:5433/molecules_test_db"
engine = create_engine(TEST_DATABASE_URL)
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Создаем тестовый клиент
client = TestClient(app)

@pytest.fixture(scope="module")
def test_db():
    models.Base.metadata.create_all(bind=engine)
    yield
    models.Base.metadata.drop_all(bind=engine)

@pytest.fixture(scope="function")
def db_session():
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.rollback()
        db.close()

@pytest.fixture(scope="function", autouse=True)
def reset_db(db_session):
    db_session.query(models.Molecule).delete()  
    db_session.add_all([
        models.Molecule(id=1, smiles="c1ccccc1"),
        models.Molecule(id=2, smiles="CCO"),
        models.Molecule(id=3, smiles="CC(=O)Oc1ccccc1C(=O)O"),
    ])
    db_session.commit()

# tests
def test_get_server():
    """Тестируем эндпоинт GET /"""
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()

def test_add_molecule():
    """Тестируем добавление молекулы"""
    molecule = {"id": 4, "smiles": "C"}
    response = client.post("/add", json=molecule)
    assert response.status_code == 200
    assert response.json()["message"] == "Molecule with id 4 added successfully."

def test_add_molecule_duplicate():
    """Тестируем добавление дублирующей молекулы"""
    molecule = {"id": 1, "smiles": "C"}
    response = client.post("/add", json=molecule)
    assert response.status_code == 400
    assert response.json()["detail"] == "Molecule already exists."

def test_add_molecule_invalid_smiles():
    """Тестируем добавление молекулы с неправильным SMILES"""
    molecule = {"id": 1000, "smiles": "invalid_smiles"}
    response = client.post("/add", json=molecule)
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid SMILES."

def test_get_molecule():
    """Тестируем получение молекулы по ID"""
    response = client.get("/search/1")
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "c1ccccc1"}

def test_get_molecule_not_found():
    """Тестируем получение несуществующей молекулы"""
    response = client.get("/search/999")
    assert response.status_code == 404
    assert response.json()["detail"] == "Molecule not found."

def test_list_molecules():
    """Тестируем список молекул"""
    response = client.get("/list")
    assert response.status_code == 200
    molecules = response.json()["molecules"]
    assert len(molecules) == 3
    assert {"id": 1, "smiles": "c1ccccc1"} in molecules

def test_update_molecule():
    """Тестируем обновление молекулы"""
    molecule = {"smiles": "CC"}
    response = client.put("/update/1", json=molecule)
    assert response.status_code == 200
    assert "Successfully updated molecule" in response.json()["message"]

def test_update_molecule_invalid_id():
    """Тестируем обновление молекулы с неверным ID"""
    molecule = {"smiles": "CC"}
    response = client.put("/update/999", json=molecule)
    assert response.status_code == 404
    assert response.json()["detail"] == "Molecule not found."

def test_update_molecule_invalid_smiles():
    """Тестируем обновление молекулы с неверным SMILES"""
    molecule = {"smiles": "invalid_smiles"}
    response = client.put("/update/1", json=molecule)
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid SMILES."

def test_delete_molecule():
    """Тестируем удаление молекулы"""
    response = client.delete("/delete/1")
    assert response.status_code == 200
    assert "Successfully deleted molecule" in response.json()["message"]

def test_delete_molecule_invalid_id():
    """Тестируем удаление молекулы с неверным ID"""
    response = client.delete("/delete/999")
    assert response.status_code == 404
    assert response.json()["detail"] == "Molecule not found."

def test_substructure_search():
    """Тестируем поиск молекулы по субструктуре"""
    response = client.get("/substructure_search/CC")
    assert response.status_code == 200
    matches = response.json()["matched_molecules"]
    assert len(matches) > 0
    for match in matches:
        assert "id" in match and "smiles" in match

def test_substructure_search_invalid_smiles():
    """Тестируем поиск молекулы по неверной субструктуре SMILES"""
    response = client.get("/substructure_search/invalid_smiles")
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid substructure SMILES."
