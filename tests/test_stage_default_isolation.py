"""Editing one block must not change other blocks or future defaults."""
import copy
import unittest

from biomolexplorer.catalog import new_stage, operation_fields


class StageDefaultIsolationTests(unittest.TestCase):
    def test_docking_box_is_independent_for_each_block(self):
        for operation in ('redocking','docking_vina'):
            with self.subTest(operation=operation):
                first, second = new_stage(operation), new_stage(operation)
                expected = copy.deepcopy(second['parameters']['sizeof_box'])
                first['parameters']['sizeof_box'][0] = 35
                self.assertEqual(second['parameters']['sizeof_box'],expected)
                third = new_stage(operation)
                self.assertEqual(third['parameters']['sizeof_box'],expected)
                default = next(f['default'] for f in operation_fields(operation)
                    if f['name']=='sizeof_box')
                self.assertEqual(default,expected)

    def test_uploaded_files_are_independent_for_each_block(self):
        first, second = new_stage('import_results'), new_stage('import_results')
        first['parameters']['asset_ids'].append('uploaded-file')
        self.assertEqual(second['parameters']['asset_ids'],[])
        self.assertEqual(new_stage('import_results')['parameters']['asset_ids'],[])


if __name__ == '__main__':
    unittest.main()
