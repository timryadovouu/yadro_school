from sqlalchemy import Column, String, Integer, Text
from database import Base
import logging 

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Molecule(Base):
    __tablename__ = "molecules"

    id = Column(String, primary_key=True, index=True)
    smiles = Column(String, nullable=False)

    def __repr__(self):
        logger.info("call repr function from Molecule class")
        return f"<Molecule(id={self.id}, smiles={self.smiles})>"
    
class TaskResult(Base):
    __tablename__ = "task_results"

    id = Column(Integer, primary_key=True, index=True)
    task_id = Column(String, unique=True, index=True)
    status = Column(String, index=True)
    result = Column(Text, nullable=True)
