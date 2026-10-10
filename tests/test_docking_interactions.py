"""Geometric classifications retain the actual docked 3D endpoints."""
import math
import unittest
from biomolexplorer.docking_interactions import interaction_diagram, ligand_molecule


def atom(serial,name,resn,x,y,z,element):
    return f'ATOM  {serial:5d} {name:<4} {resn:>3} A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n'


class InteractionTests(unittest.TestCase):
    def molecule(self,smiles,coordinates):
        from rdkit import Chem
        mol=Chem.MolFromSmiles(smiles)
        conf=Chem.Conformer(mol.GetNumAtoms())
        for i,xyz in enumerate(coordinates):conf.SetAtomPosition(i,xyz)
        mol.AddConformer(conf)
        return {'format':'sdf','data':Chem.MolToMolBlock(mol)}

    def test_parallel_and_perpendicular_aromatic_rings(self):
        coordinates=[(1.4*math.cos(i*math.pi/3),1.4*math.sin(i*math.pi/3),0) for i in range(6)]
        pose=self.molecule('c1ccccc1',coordinates)
        names=('CG','CD1','CE1','CZ','CE2','CD2')
        for perpendicular,kind in ((False,'pi_parallel'),(True,'pi_t')):
            points=[(x,0,y+3.5) if perpendicular else (x,y,3.5) for x,y,_ in coordinates]
            receptor={'format':'pdb','data':''.join(atom(i+1,n,'PHE',*p,'C') for i,(n,p) in enumerate(zip(names,points)))}
            result=interaction_diagram(receptor,pose)
            self.assertTrue(result['available'])
            row=next(r for r in result['interactions'] if r['kind']==kind)
            self.assertAlmostEqual(row['distance'],3.5,places=2)
            self.assertAlmostEqual(row['ligand_position'][2],0)
            self.assertAlmostEqual(row['receptor_position'][2],3.5,places=2)

    def test_hydrogen_bond_requires_explicit_h_and_correct_geometry(self):
        pose=self.molecule('CO',[(0,0,0),(1.4,0,0)])
        for hydrogen,expected in ((None,False),((3.2,0,0),True),((4.2,1,0),False)):
            text=atom(1,'OG','SER',4.2,0,0,'O')
            if hydrogen:text+=atom(2,'HG','SER',*hydrogen,'H')
            result=interaction_diagram({'format':'pdb','data':text},pose)
            bonds=[r for r in result['interactions'] if r['kind']=='hydrogen_bond']
            self.assertEqual(bool(bonds),expected)
            if expected:self.assertEqual(bonds[0]['receptor_position'],[4.2,0,0])

    def test_vdw_rejects_overlap_and_distant_atoms(self):
        pose=self.molecule('C',[(0,0,0)])
        for distance,expected in ((3.4,True),(.8,False),(5,False)):
            receptor={'format':'pdb','data':atom(1,'CB','ALA',distance,0,0,'C')}
            result=interaction_diagram(receptor,pose)
            self.assertEqual(any(r['kind']=='van_der_waals' for r in result['interactions']),expected)

    def test_vina_topology_accepts_present_hydrogens_without_moving_atoms(self):
        text=atom(1,'C1','LIG',0,0,0,'C')+atom(2,'H1','LIG',1.1,0,0,'H')
        # PDBQT atom type resides at the end of the record.
        text='\n'.join(line[:70]+' 0.000 '+('HD' if 'H1' in line else 'C') for line in text.splitlines())
        mol=ligand_molecule({'format':'pdbqt','data':text},'C')
        self.assertIsNotNone(mol)
        self.assertEqual(mol.GetNumAtoms(),1)
        self.assertEqual(list(mol.GetConformer().GetAtomPosition(0)),[0,0,0])

    def test_ambiguous_topology_disables_classification(self):
        result=interaction_diagram({'format':'pdb','data':atom(1,'CB','ALA',3,0,0,'C')},
            {'format':'pdbqt','data':atom(1,'C1','LIG',0,0,0,'C')},'CCO')
        self.assertFalse(result['available'])
