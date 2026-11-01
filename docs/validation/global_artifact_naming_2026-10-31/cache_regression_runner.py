import unittest,json
from pathlib import Path
modules=['test_stage_cache','test_restart_reuse','test_analysis_origins','test_docking_incremental','test_redocking_incremental','test_dock6_incremental']
# Exact failing test IDs reproduced against untouched HEAD in a /tmp copy.
import re
text=Path(__file__).with_name('baseline-regression.log').read_text()
excluded=set(re.findall(r'^(?:FAIL|ERROR): (test_\w+) ',text,re.M))
def flat(suite):
 for test in suite:
  if isinstance(test,unittest.TestSuite):yield from flat(test)
  else:yield test
all_tests=list(flat(unittest.defaultTestLoader.loadTestsFromNames(modules)))
selected=[t for t in all_tests if t._testMethodName not in excluded]
print('Baseline exclusions:',json.dumps(sorted(excluded)),flush=True)
result=unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(selected));raise SystemExit(not result.wasSuccessful())
