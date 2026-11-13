"""Generic PTM/structural contracts independent of any target identifier."""
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from biomolexplorer.md_residue_support import bundle,plan,pdb_residues
from biomolexplorer.md_torsion_support import read_torsion_parameters
from biomolexplorer.structure_sanitation import sanitize,receptor_script


def atom(serial,name,residue,number,x,y=0,z=0,alt='',occupancy=1,record='HETATM',element=None):
    element=element or name[0]
    return f'{record:6s}{serial:5d} {name:>4s}{alt or " "}{residue:3s} A{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{20:6.2f}          {element:>2s}\n'


class GenericResidueTests(unittest.TestCase):
    def test_all_published_models_have_integrity_graphs_and_charge(self):
        data=bundle();self.assertEqual(len(data['residues']),32)
        for code,model in data['residues'].items():
            with self.subTest(code=code):
                self.assertTrue({'N','CA','C','O'}<=set(model['atoms']))
                self.assertAlmostEqual(sum(a['charge'] for a in model['atoms'].values()),model['net_charge'],places=3)
                self.assertTrue(all(set(edge)<=set(model['atoms']) for edge in model['bonds']))
        self.assertEqual(data['residues']['CSS']['net_charge'],0)
        for code in ('TPO','SEP','PTR'):self.assertEqual(data['residues'][code]['net_charge'],-2)

    def test_altloc_is_residue_consistent_and_source_immutable(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'input.pdb';target=Path(temp)/'output.pdb'
            text=atom(1,'N','ALA',1,0,record='ATOM')+atom(2,'CA','ALA',1,1,alt='A',occupancy=.2,record='ATOM')+atom(3,'CA','ALA',1,2,alt='B',occupancy=.8,record='ATOM')
            source.write_text(text);report=sanitize(source,target)
            self.assertEqual(source.read_text(),text)
            residues,serials,_=pdb_residues(target)
            self.assertEqual(residues['A   1 ']['atoms']['CA']['xyz'][0],2)
            self.assertEqual(set(serials),{1,3});self.assertEqual(report['alternate_conformations'][0]['selected'],'B')

    def test_cofactor_and_covalent_ptm_are_not_deleted_as_additives(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'input.pdb';target=Path(temp)/'output.pdb'
            source.write_text(atom(1,'C','ALA',1,0,record='ATOM')+atom(2,'N','TPO',2,1.3)+atom(3,'C1','FAD',3,15)+atom(4,'O','HOH',4,20)+atom(5,'N','GLY',5,30)+atom(6,'C1','XYZ',6,40))
            report=sanitize(source,target);roles={r['name']:r for r in report['residue_roles']}
            self.assertEqual(roles['TPO']['role'],'covalently_connected')
            for name in ('TPO','FAD','XYZ'):self.assertEqual(roles[name]['action'],'retained')
            for name in ('HOH','GLY'):self.assertEqual(roles[name]['action'],'removed')

    def test_explicit_cofactor_overrides_additive_removal(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'in.pdb';source.write_text(atom(1,'O','GOL',1,1))
            report=sanitize(source,Path(temp)/'out.pdb',keep=['GOL'])
            self.assertEqual(report['residue_roles'][0]['action'],'retained')

    def test_redocking_preserves_cofactors_but_excludes_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'in.pdb';source.write_text(atom(1,'CA','ALA',1,0,record='ATOM')+atom(2,'C1','FAD',2,10)+atom(3,'C1','LIG',3,20))
            text='open '+str(source)+'\nselect protein\nselect invert\ndelete selected\naddh\n'
            updated=receptor_script(text,Path(temp)/'prep.com',['9XYZ','LIG',3,'A'])
            self.assertIn('| :FAD',updated)
            sanitized=Path(temp)/'prep.sanitized.pdb';self.assertNotIn('LIG',sanitized.read_text());self.assertIn('FAD',sanitized.read_text())

    def test_multi_model_duplicate_and_nonfinite_coordinates_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'in.pdb';target=Path(temp)/'out.pdb'
            for text in ('MODEL        1\nMODEL        2\n',atom(1,'CA','ALA',1,0)+atom(2,'CA','ALA',1,1),atom(1,'CA','ALA',1,float('nan'))):
                with self.subTest(text=text):
                    source.write_text(text)
                    with self.assertRaises(ValueError):sanitize(source,target)

    def test_native_profile_does_not_silently_rewrite_unknown_residue(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'input.pdb';path.write_text(atom(1,'SE','MSE',1,0,element='SE'))
            report=plan(path);self.assertEqual(report['profile'],'native_ff94');self.assertIn('MSE',path.read_text())

    def test_torsion_scan_checks_arbitrary_phases_and_periodicities(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'top.prmtop'
            def topology(numbers,phases):
                path.write_text('%FLAG DIHEDRAL_PERIODICITY\n%FORMAT(5E16.8)\n'+''.join(f'{v:16.8E}' for v in numbers)+'\n%FLAG DIHEDRAL_PHASE\n%FORMAT(5E16.8)\n'+''.join(f'{v:16.8E}' for v in phases)+'\n')
            topology([1,2,3,4,6],[0,3.141594,0,math.pi,0]);self.assertEqual(read_torsion_parameters(path),[])
            topology([5,-7,10],[.5,1.3,math.pi]);self.assertEqual(len(read_torsion_parameters(path)),3)
            topology([1.5],[0])
            with self.assertRaises(ValueError):read_torsion_parameters(path)

    def test_unknown_covalent_component_never_receives_detached_gaff_charges(self):
        from biomolexplorer.md_cofactor_support import parameterize
        with tempfile.TemporaryDirectory() as temp:
            residues={'A   1 ':dict(name='MSE',atoms={name:dict(element=element,xyz=(0,0,0)) for name,element in [('N','N'),('CA','C'),('C','C'),('SE','SE')]})}
            with patch('biomolexplorer.md_cofactor_support.native_residues',return_value=set()),patch('biomolexplorer.preparation_tools.run_amber_tool',side_effect=AssertionError('Unsafe chemistry')):
                with self.assertRaisesRegex(ValueError,'bonded protein profile'):parameterize(temp,temp,residues,{})

    def test_valid_dynamic_checkpoint_requires_no_parameter_tool(self):
        from biomolexplorer.md_cofactor_support import parameterize
        from biomolexplorer.stage_cache import file_digest
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'XYZ.mol2').write_text('verified');(root/'XYZ.frcmod').write_text('verified')
            report={'dynamic_cofactors':[dict(component='XYZ',files={name:file_digest(root/name) for name in ('XYZ.mol2','XYZ.frcmod')})]}
            with patch('biomolexplorer.preparation_tools.resolve_amber_tool',side_effect=AssertionError('Cache should skip tools')):
                self.assertEqual(len(parameterize(root,root,{},report)),2)


class NativePeriodicDerivativeTests(unittest.TestCase):
    def test_general_torsion_energy_and_gradient_against_independent_dihedral(self):
        import ctypes
        import random
        import shutil
        import subprocess
        from rdkit import Chem
        from rdkit.Chem import rdMolTransforms
        from biomolexplorer.md_torsion_support import HEADER
        compiler=shutil.which('gcc')
        if not compiler:self.skipTest('Native C compiler unavailable')
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);shutil.copy2(HEADER,root/HEADER.name)
            source=root/'test.c';source.write_text('#include "periodic_torsion.h"\ndouble evaluate(double *x,double n,double k,double phase,double *g) { return bme_periodic_torsion(x,n,k,phase,g); }\n')
            library=root/'test.so';subprocess.run([compiler,'-shared','-fPIC','-O2',str(source),'-o',str(library),'-lm'],check=True,capture_output=True)
            evaluate=ctypes.CDLL(str(library)).evaluate;array=ctypes.c_double*12
            evaluate.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)];evaluate.restype=ctypes.c_double
            randomizer=random.Random(3416);mol=Chem.MolFromSmiles('CCCC');conformer=Chem.Conformer(4)
            for sample in range(100):
                coords=[randomizer.uniform(-3,3) for _ in range(12)]
                n=randomizer.randint(1,12);phase=randomizer.uniform(-math.pi,math.pi);k=randomizer.uniform(.1,3)
                for index in range(4):conformer.SetAtomPosition(index,coords[index*3:index*3+3])
                phi=rdMolTransforms.GetDihedralRad(conformer,0,1,2,3);gradient=array()
                value=evaluate(array(*coords),n,k,phase,gradient)
                self.assertAlmostEqual(value,k*(1+math.cos(n*phi-phase)),places=10)
                for index in range(12):
                    plus=list(coords);minus=list(coords);plus[index]+=1e-6;minus[index]-=1e-6
                    finite=(evaluate(array(*plus),n,k,phase,array())-evaluate(array(*minus),n,k,phase,array()))/2e-6
                    self.assertAlmostEqual(gradient[index],finite,delta=2e-6)


class DeferredRedockingSanitationTests(unittest.TestCase):
    def test_receptor_is_sanitized_after_complex_becomes_available(self):
        from biomolexplorer.chimera_preparation import receptor_cleanup
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'late.complex.pdb';script=root/'prepare_receptor_X_A.com'
            text='open '+str(source)+'\nselect protein\nselect invert\ndelete selected\naddh\naddcharge all spec #0 chargeModel 14sb method gas\nwrite format mol2 #0 '+str(root/'rec.mol2')+'\n'
            self.assertEqual(receptor_script(text,script),text)
            script.write_text(text)
            script.with_suffix('.receptor_selection.json').write_text(json.dumps({'record':['9XYZ','LIG',3,'A']}))
            source.write_text(atom(1,'CA','ALA',1,0,record='ATOM')+atom(2,'C1','FAD',2,10)+atom(3,'C1','LIG',3,20)+atom(4,'O','HOH',4,30))
            receptor_cleanup(script)
            sanitized=script.with_suffix('.sanitized.pdb')
            self.assertTrue(sanitized.is_file());self.assertIn('FAD',sanitized.read_text());self.assertNotIn('LIG',sanitized.read_text());self.assertNotIn('HOH',sanitized.read_text())
            self.assertIn('open '+str(sanitized),script.read_text());self.assertIn('| :FAD',script.read_text())


class EmptyPreparationLogCacheTests(unittest.TestCase):
    def test_completed_pair_with_reused_preparation_empty_logs_is_cached(self):
        from biomolexplorer.md_cache import PairCache
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp);directory=output/'pairs'/'one';directory.mkdir(parents=True)
            (directory/'analysis.json').write_text('{"validated": true}')
            (directory/'amberize_ligand.lig.1.out').touch()
            record=dict(receptor='target',ligand='compound',status='succeeded',amber_score=-1)
            cache=PairCache(output);cache.remember('completed',directory,record)
            restored=cache.restore('completed',output/'repeat')
            self.assertEqual(restored['status'],'succeeded')
            self.assertTrue((output/'repeat/analysis.json').is_file())
            self.assertFalse((output/'repeat/amberize_ligand.lig.1.out').exists())

    def test_empty_scientific_topology_is_still_rejected(self):
        from biomolexplorer.md_cache import PairCache
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp);directory=output/'pairs'/'one';directory.mkdir(parents=True)
            (directory/'rec.prmtop').touch()
            with self.assertRaises(ValueError):
                PairCache(output).remember('invalid',directory,dict(receptor='target',ligand='compound'))
            self.assertFalse((output/'md_cache.json').exists())


if __name__=='__main__':unittest.main()
