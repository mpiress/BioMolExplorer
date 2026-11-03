"""Offline packaging smoke tests; run with the scientific environment's Python.

These tests reuse installed dependencies, not a fresh Conda solve.
"""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import tomllib
import unittest
import venv

ROOT = Path(__file__).resolve().parents[1]


class InstallationGuideTests(unittest.TestCase):
    def test_bilingual_commands_and_shell_syntax(self):
        for path in (ROOT / 'docs/installation.md', ROOT / 'docs/en/installation.md'):
            text = path.read_text()
            self.assertIn("python -m pip install -e '.[ui]'", text)
            self.assertLess(text.index("python -m pip install -e '.[ui]'"),
                            text.index('biomolexplorer-ui --help'))
            self.assertIn('in chimera dock6 sphgen showbox grid obabel vina;', text)
            for block in re.findall(r'```bash\n(.*?)```', text, re.S):
                subprocess.run(['bash', '-n'], input=block, text=True, check=True)

    def test_external_tools_are_not_pip_dependencies(self):
        config = tomllib.loads((ROOT / 'pyproject.toml').read_text())
        self.assertEqual(config['project']['dependencies'], [])
        self.assertIn('flet[desktop]==1.0.3', config['project']['optional-dependencies']['ui'])
        self.assertNotIn('admet-ai==', (ROOT / 'requirements.yml').read_text())

    def test_editable_install_registers_working_commands_outside_checkout(self):
        with tempfile.TemporaryDirectory(prefix='biomol-install-') as directory:
            target = Path(directory)
            venv.create(target / 'env', with_pip=True, system_site_packages=True)
            python = target / 'env/bin/python'
            env = dict(os.environ, PIP_DISABLE_PIP_VERSION_CHECK='1')
            env.pop('PYTHONPATH', None)
            subprocess.run([str(python), '-m', 'pip', 'install', '--no-deps',
                            '--no-build-isolation', '-e', str(ROOT) + '[ui]'],
                           env=env, cwd=target, check=True, capture_output=True, text=True)
            for name in ('biomolexplorer', 'biomolexplorer-ui'):
                result = subprocess.run([str(target / 'env/bin' / name), '--help'],
                                        env=env, cwd=target, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn('usage:', result.stdout)
            result = subprocess.run([str(target / 'env/bin/biomolexplorer-ui'), '--help'],
                                    env=env, cwd=target, capture_output=True, text=True)
            for flag in ('--web', '--language', '--dock6-path', '--no-browser'):
                self.assertIn(flag, result.stdout)


if __name__ == '__main__':
    unittest.main()
