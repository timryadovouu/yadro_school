from celery import Celery 
import tasks

celery = Celery( 
    'tasks', 
    broker='redis://redis:6379/0', 
    backend='redis://redis:6379/0',
) 

celery.conf.update(task_track_started=True)
celery.autodiscover_tasks(['tasks'])