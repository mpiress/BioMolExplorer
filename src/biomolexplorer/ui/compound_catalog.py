"""Per-dialog, on-demand references; scientific validation stays off the UI loop."""
import asyncio
from collections import OrderedDict
from pathlib import Path


class CompoundCatalog:
    def __init__(self, ui, paths, valid):
        self.ui, self.paths, self.valid = ui, paths, valid
        self.references = {}
        self.lock = asyncio.Lock()

    async def load(self, origin):
        if not origin or origin not in self.paths:
            return
        async with self.lock:
            if not self.valid():
                return
            records = {}
            for path in self.paths[origin]:
                # Cache mutation belongs to the event loop. A worker finishing
                # after logout must never repopulate the session's private cache.
                key = await self.ui.call(self._signature, path)
                if not self.valid():
                    return
                cache = getattr(self.ui, '_compound_reference_cache', None)
                if cache is None:
                    cache = self.ui._compound_reference_cache = OrderedDict()
                if key in cache:
                    cache.move_to_end(key)
                    values = cache[key]
                else:
                    from biomolexplorer.pubchem_retrieval import reference_choices
                    values = await self.ui.call(reference_choices, [path])
                    if not self.valid():
                        return
                    # A file changed during reading cannot enter the cache.
                    current = await self.ui.call(self._signature, path)
                    if not self.valid():
                        return
                    if key == current and len(values) <= 20000:
                        for previous in list(cache):
                            if previous[0] == str(path):
                                del cache[previous]
                        cache[key] = values
                        while len(cache) > 32 or sum(len(v) for v in cache.values()) > 20000:
                            cache.popitem(last=False)
                records[str(path)] = values
            self.references[origin] = records

    @staticmethod
    def _signature(path):
        path = Path(path)
        stat = path.stat()
        return (str(path), stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size)
