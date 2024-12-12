from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import os
import logging

# DATABASE_URL =os.getenv("DATABASE_URL")
DATABASE_URL = os.getenv("DATABASE_URL", f"postgresql://${os.getenv("POSTGRES_USER")}:${os.getenv("POSTGRES_PASSWORD")}@db:5432/${os.getenv("POSTGRES_DB")}")

engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_db():
    logger.info("call get_db function")
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()