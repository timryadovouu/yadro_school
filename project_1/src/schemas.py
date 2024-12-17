from typing import Optional
from pydantic import BaseModel

class TaskResultBase(BaseModel):
    task_id: str
    status: str
    result: Optional[str] = None

    class Config:
        orm_mode = True