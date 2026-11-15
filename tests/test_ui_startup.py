import unittest

from biomolexplorer.ui.app import build_parser, resolve_web_startup


class UIStartupTests(unittest.TestCase):
    def startup(self, *arguments):
        return resolve_web_startup(build_parser().parse_args(arguments))

    def test_web_remains_local_only_by_default(self):
        self.assertEqual(self.startup('--web'), (True, '127.0.0.1'))

    def test_web_server_binds_all_network_interfaces(self):
        self.assertEqual(self.startup('--web-server'), (True, '0.0.0.0'))

    def test_explicit_host_overrides_mode_default(self):
        self.assertEqual(self.startup('--web-server', '--host', '192.168.1.20'),
                         (True, '192.168.1.20'))
        self.assertEqual(self.startup('--web', '--host', '0.0.0.0'),
                         (True, '0.0.0.0'))

    def test_no_browser_keeps_legacy_implicit_web_mode(self):
        self.assertEqual(self.startup('--no-browser'), (True, '127.0.0.1'))


if __name__ == '__main__':
    unittest.main()
