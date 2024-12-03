from sqlalchemy import Column, Integer, String
import database

class Molecule(database.Base):
    __tablename__ = "molecules"

    id = Column(String, primary_key=True, index=True)
    smiles = Column(String, nullable=False)

    def __repr__(self):
        return f"<Molecule(id={self.id}, smiles={self.smiles})>"
