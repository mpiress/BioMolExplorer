"""Mathematical invariants of the compound-level docking consensus."""
import json
import math
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from biomolexplorer.consensus_model import consensus_rows


def record(code, target, score, smiles='CCO'):
    return dict(molecule_chembl_id=code, receptor_id=target, score=score, canonical_smiles=smiles)


class ConsensusModelTests(unittest.TestCase):
    def test_single_target_closed_form_and_positive_energies(self):
        vina = [record('A', 'T', -2), record('B', 'T', 1), record('C', 'T', 100)]
        dock6 = [record('A', 'T', -20), record('B', 'T', -10), record('C', 'T', 5)]
        rows = consensus_rows(vina, dock6)
        for row, expected in zip(rows, (1., math.exp(-5), math.exp(-10))):
            self.assertAlmostEqual(row['vina'], expected)
            self.assertAlmostEqual(row['dock6'], expected)
            self.assertAlmostEqual(row['normalized_score'], expected)
        self.assertEqual(rows[0]['vina_raw'], -2)

    def test_monotone_scale_and_offset_invariance(self):
        vina = [record(str(i), 'T', i) for i in range(10)]
        original = consensus_rows(vina, vina)
        transformed = consensus_rows([dict(r, score=7*r['score']+200) for r in vina], vina)
        self.assertEqual([r['consensus'] for r in original], [r['consensus'] for r in transformed])

    def test_1000_compounds_three_targets_one_row_each(self):
        data = [record(f'M{i:04}', target, i) for target in ('T1', 'T2', 'T3') for i in range(1000)]
        rows = consensus_rows(data, data)
        self.assertEqual(len(rows), 1000)
        self.assertEqual(len({r['molecule_chembl_id'] for r in rows}), 1000)
        self.assertTrue(all(r['paired_target_count'] == 3 for r in rows))
        self.assertTrue(all(len(json.loads(r['target_scores'])) == 3 for r in rows))

    def test_breadth_penalizes_one_target_specialist(self):
        data = [record(code, target, value) for target, values in (
            ('T1', (0, 1, 2)), ('T2', (2, 0, 1)), ('T3', (2, 0, 1)))
                for code, value in zip(('specialist', 'broad', 'other'), values)]
        rows = consensus_rows(data, data)
        self.assertEqual(rows[0]['molecule_chembl_id'], 'broad')

    def test_missing_target_is_penalized_and_not_silently_dropped(self):
        data = [record('A', 'T1', -2), record('B', 'T1', -1), record('B', 'T2', -1)]
        rows = {r['molecule_chembl_id']: r for r in consensus_rows(data, data)}
        self.assertEqual(rows['A']['vina_target_coverage'], .5)
        self.assertAlmostEqual(rows['A']['vina'], .5*math.exp(-5))
        self.assertIsNone(json.loads(rows['A']['target_scores'])[1]['vina_raw'])

    def test_duplicate_poses_do_not_change_distribution_or_result(self):
        data = [record('A', 'T', -2), record('B', 'T', -1)]
        self.assertEqual(consensus_rows(data, data), consensus_rows(data*3+[record('A', 'T', 100)], data*4))

    def test_ties_and_singleton_are_finite_and_deterministic(self):
        data = [record(code, 'T', -1) for code in ('C', 'A', 'B')]
        rows = consensus_rows(data, data)
        self.assertEqual([r['molecule_chembl_id'] for r in rows], ['A', 'B', 'C'])
        self.assertTrue(all(abs(r['consensus']-math.exp(-5)) < 1e-12 for r in rows))
        self.assertEqual(consensus_rows(data[:1], data[:1])[0]['consensus'], 1.)

    def test_input_order_and_target_order_invariance(self):
        data = [record(str(i), target, i) for target in ('T2', 'T1') for i in range(10)]
        expected = consensus_rows(data, data)
        random.Random(11).shuffle(data)
        self.assertEqual(consensus_rows(data, data), expected)

    def test_stereoisomers_never_merge_under_same_identifier(self):
        a = record('M', 'T1', -1, 'N[C@@H](C)C(=O)O')
        b = record('M', 'T2', -1, 'N[C@H](C)C(=O)O')
        with self.assertRaisesRegex(ValueError, 'estruturas diferentes'):
            consensus_rows([a, b], [a])
        b['molecule_chembl_id'] = 'other'
        self.assertEqual(len(consensus_rows([a, b], [a, b])), 2)

    def test_invalid_scores_fail_and_different_targets_never_pair(self):
        for value in (math.nan, math.inf, -math.inf):
            with self.assertRaises(ValueError):
                consensus_rows([record('A', 'T', value)], [record('A', 'T', -1)])
        self.assertEqual(consensus_rows([record('A', 'T1', -1)], [record('A', 'T2', -1)]), [])

    def test_pose_and_context_share_one_representative_target(self):
        data = [dict(record('A', t, -1), conformer_file=t+'.mol2', receptor_file=t+'.pdb') for t in ('T2', 'T1')]
        row = consensus_rows(data, data)[0]
        self.assertEqual(row['receptor_id'], 'T1')
        self.assertEqual(row['vina_pose'], 'T1.mol2')
        self.assertEqual(row['dock6_receptor_file'], 'T1.pdb')

    def test_export_has_exactly_one_thousand_rows_and_target_audit(self):
        from biomolexplorer.docking_data import write_csv, read_compounds
        from wrappers.docking import generate_consensus
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for engine in ('vina', 'dock6'):
                folder = root/engine
                folder.mkdir()
                data = [record(str(i), target, i) for target in ('T1', 'T2', 'T3') for i in range(1000)]
                write_csv(folder/'docking_results.csv', data)
            with patch('wrappers.docking.plot_scatter_comparison') as plot:
                result = generate_consensus(str(root), str(root/'out'), 'ranking',
                    base_vina_path=str(root/'vina'), base_dock6_path=str(root/'dock6'))
            plot.assert_called_once()
            exported = read_compounds(root/'out/ranking.csv')
            self.assertEqual(len(result), 1000)
            self.assertEqual(len(exported), 1000)
            self.assertTrue(all(len(json.loads(row['target_scores'])) == 3 for row in exported))
            self.assertEqual(exported[0]['consensus_model'], 'target-balanced-ecr-v1')
            self.assertEqual(float(exported[0]['normalized_score']), 1.)

    def test_consensus_identity_is_distinct_from_shared_docking_wrapper(self):
        from biomolexplorer.stage_cache import operation_digest
        self.assertNotEqual(operation_digest('consensus'), operation_digest('docking_vina'))
        self.assertNotEqual(operation_digest('consensus'), operation_digest('docking_dock6'))

    def test_new_labels_and_validation_errors_translate(self):
        from biomolexplorer.ui.localization import Translator
        translate = Translator('en')
        for source, expected in (
            ('Relevância Vina', 'Vina relevance'),
            ('Relevância DOCK6', 'DOCK6 relevance'),
            ('Consenso', 'Consensus'),
            ('Consenso: código molecular vazio.', 'Consensus: empty compound identifier.'),
            ('Consenso: score não finito para M1', 'Consensus: non-finite score for M1'),
        ):
            self.assertEqual(translate(source), expected)


if __name__ == '__main__':
    unittest.main()
