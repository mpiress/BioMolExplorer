"""Coordinate provider requests across local worker processes (Linux)."""
import time
from pathlib import Path


def wait_for_slot(path, interval=0.25):
    import fcntl
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        stream.seek(0)
        last = float(stream.read() or 0)
        now = time.monotonic()
        # After a machine reboot, monotonic timestamps may move backwards.
        delay = max(0, interval - (now - last)) if last <= now else 0
        time.sleep(delay)
        stream.seek(0)
        stream.truncate()
        stream.write(str(time.monotonic()))
        stream.flush()
