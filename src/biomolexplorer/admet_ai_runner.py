"""JSON boundary for the isolated optional predictor; logs go to stderr."""
import contextlib
import json
import sys
from biomolexplorer.admet_predictors import predict_local

if __name__=='__main__':
    smiles=json.load(sys.stdin)
    with contextlib.redirect_stdout(sys.stderr):
        records,provider=predict_local(smiles)
    json.dump(dict(records=records,provider=provider),sys.stdout,allow_nan=False)
