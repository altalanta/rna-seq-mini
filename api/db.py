import os
from sqlalchemy import create_engine, Column, String, DateTime, Text
from sqlalchemy.orm import sessionmaker, declarative_base
import datetime

# Database URL - configurable via environment variable
# Supports SQLite (default), PostgreSQL, MySQL, etc.
DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///./api_jobs.db")

# SQLite requires check_same_thread=False for FastAPI compatibility
connect_args = {"check_same_thread": False} if "sqlite" in DATABASE_URL else {}
engine = create_engine(DATABASE_URL, connect_args=connect_args)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

class Job(Base):
    __tablename__ = "jobs"

    id = Column(String, primary_key=True, index=True)
    status = Column(String, default="queued")
    created_at = Column(DateTime, default=datetime.datetime.utcnow)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    parameters = Column(Text, nullable=True)
    results_dir = Column(String, nullable=True)
    log_file = Column(String, nullable=True)
    error = Column(Text, nullable=True)

def init_db():
    Base.metadata.create_all(bind=engine)

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


