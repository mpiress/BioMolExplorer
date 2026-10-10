"""Include the installed commit in wheels, which do not retain Git metadata."""
import subprocess
import re
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py


class BuildWithRevision(build_py):
    def run(self):
        super().run()
        root = Path(__file__).parent
        revision = None
        if (root / '.git').exists():
            revision = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'],
                capture_output=True, text=True, check=True).stdout.strip()
        else:
            match = re.fullmatch(r'BioMolExplorer-([a-f0-9]{40})', root.resolve().name)
            if match: revision = match[1]
        if revision:
            output = Path(self.build_lib) / 'biomolexplorer/resources/revision.txt'
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(revision + '\n')


setup(cmdclass={'build_py': BuildWithRevision})
