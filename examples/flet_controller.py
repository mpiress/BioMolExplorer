"""UI-facing async controller. Add Flet controls without importing science here.

Use from the Flet event loop. UI code receives only serializable job snapshots;
all SQLite calls and process cancellation run outside that loop.
"""
import asyncio
from biomolexplorer.jobs import JobManager, TERMINAL


class WorkflowController:
    def __init__(self, config):
        self.jobs = JobManager(config)

    async def submit(self, operation, parameters):
        return await asyncio.to_thread(self.jobs.submit, operation, parameters)

    async def get(self, job_id):
        return await asyncio.to_thread(self.jobs.get, job_id)

    async def cancel(self, job_id):
        return await asyncio.to_thread(self.jobs.cancel, job_id)

    async def updates(self, job_id, interval=0.5):
        while True:
            job = await self.get(job_id)
            yield job
            if job['status'] in TERMINAL:
                return
            await asyncio.sleep(interval)

    async def close(self):
        await asyncio.to_thread(self.jobs.close)
