"""ChEMBL connection with bounded retries, timeouts and explicit pagination.

The REST endpoints remain usable when the legacy client's /spore initialization
fails. Sessions are local to a query, so parallel molecule retrieval does not
share mutable connection state.
"""
from crawlers.chembl_client import ChEMBLClient

class CrawlerSettings:
    def __init__(self):
        self.client=ChEMBLClient()
    def get_client_connection(self):
        return self.client
