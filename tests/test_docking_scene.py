"""Coordinate fidelity, pose ranking, contact geometry and capability authorization."""
import json
import tempfile
import unittest
from pathlib import Path
from biomolexplorer.docking_scene import best_pose,atoms,contacts,scene_payload
from biomolexplorer.pdb_view import StructureViewers,PREFIX
from biomolexplorer.workspace import WorkspaceStore, AccessDenied
from biomolexplorer.catalog import new_stage


def atom(serial,name,resn,chain,resi,x,y=0,z=0,element='C',record='ATOM',icode=''):
    return f'{record:<6}{serial:5d} {name:<4} {resn:>3} {chain}{resi:4d}{icode:1}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n'


class ContactGeometryTests(unittest.TestCase):
    def test_heavy_atoms_cutoff_chains_and_insertion_codes(self):
        receptor=atom(1,'CA','ALA','A',1,0)+atom(2,'H','GLY','B',2,2,element='H')+atom(3,'CA','TYR','B',3,4,icode='A')+atom(4,'O','HOH','A',4,1,element='O',record='HETATM')
        ligand=atom(5,'C1','LIG','L',1,0,record='HETATM')
        rows=contacts(atoms(receptor,'pdb'),atoms(ligand,'pdb'),4)
        self.assertEqual([(r['resn'],r['chain'],r['resi'],r['icode'],r['distance']) for r in rows],
            [('ALA','A',1,'',0.),('TYR','B',3,'A',4.)])
        self.assertEqual(len(contacts(atoms(receptor,'pdb'),atoms(ligand,'pdb'),3.9)),1)
        for value in (1,9,float('nan')):
            with self.assertRaises(ValueError):contacts([],[],value)

    def test_vina_chooses_lowest_score_without_transforming_coordinates(self):
        text='MODEL 1\nREMARK VINA RESULT: -2 0 0\n'+atom(1,'C1','LIG','A',1,8)+'ENDMDL\nMODEL 2\nREMARK VINA RESULT: -8 0 0\n'+atom(1,'C1','LIG','A',1,3)+'ENDMDL\n'
        selected,info=best_pose(text,'pdbqt')
        self.assertEqual(info['model'],2);self.assertEqual(info['score'],-8)
        self.assertIn('   3.000',selected);self.assertNotIn('   8.000',selected)

    def test_dock6_chooses_lowest_score_with_its_mol2(self):
        def molecule(score,x):return f'########## Name: LIG\n########## Grid_Score: {score}\n@<TRIPOS>MOLECULE\nLIG\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 {x} 0 0 C.3 1 LIG1 0\n@<TRIPOS>BOND\n'
        selected,info=best_pose(molecule(-1,8)+molecule(-7,3),'mol2')
        self.assertEqual(info['model'],2);self.assertEqual(info['score'],-7)
        self.assertEqual(atoms(selected,'mol2')[0]['xyz'],[3.,0.,0.])


class SceneCapabilityTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.store=WorkspaceStore(self.tmp.name);self.token=self.store.register('Owner','scene@example.org','scene-test-password')
        self.pid=self.store.create_project(self.token,'Scenes')['id'];self.root=self.store.project_dir(self.pid)
        self.rid='scene-run';self.stage=new_stage('redocking');self.sid=self.stage['id']
        with self.store.connect() as db:
            db.execute('INSERT INTO runs(id,project_id,user_id,status,stages,created,updated) VALUES(?,?,?,?,?,?,?)',
                (self.rid,self.pid,self.store.user(self.token)['id'],'succeeded','[]',1,1))
        self.viewers=StructureViewers(self.store);self.addCleanup(self.viewers.close)

    def save(self,operation,files,**extras):
        stage=dict(id=self.sid,operation=operation,status='succeeded',configuration=self.stage,artifacts=[str(p) for p in files],**extras)
        with self.store.connect() as db:db.execute('UPDATE runs SET stages=? WHERE id=?',(json.dumps([stage]),self.rid))

    def redocking(self):
        folder=self.root/'artifact/structures/Target';prepared=folder/'Prepared';prepared.mkdir(parents=True)
        out=self.root/'artifact/Target';out.mkdir()
        pdb=folder/'1ABC.pdb';pdb.write_text(atom(1,'CA','ALA','A',2,0)+atom(2,'C1','LIG','A',1,2,record='HETATM'))
        pose=out/'1ABC_LIG_1A.lig.pdbqt';pose.write_text('MODEL 1\nREMARK VINA RESULT: -2 0 0\n'+atom(2,'C1','LIG','A',1,3)+'ENDMDL\n')
        metadata=folder/'pdb_codes.csv';metadata.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RMSD\n1ABC,LIG,1,A,1\n')
        self.save('redocking',[pdb,pose,metadata]);return pdb,pose

    def test_overlay_uses_crystal_ligand_and_revokes_all_layers(self):
        from biomolexplorer.redocking_results import RedockingResults
        pdb,pose=self.redocking();service=RedockingResults(self.store)
        key=service.simulations(self.token,self.pid,self.rid,self.sid)[0]['id']
        spec=service.scene(self.token,self.pid,self.rid,self.sid,key)
        payload=scene_payload(self.store,self.token,self.pid,spec)
        self.assertEqual(payload['reference_kind'],'crystal')
        self.assertEqual(payload['contacts'][0]['distance'],3.)
        self.assertEqual(payload['reference_contacts'][0]['distance'],2.)
        self.assertNotIn('LIG',payload['layers'][0]['data'])
        ticket=self.viewers.issue_result(self.token,self.pid,self.rid,self.sid,'redocking',key,'en')
        status,mime,data,_=self.viewers.response(PREFIX+'/'+ticket+'/scene');self.assertEqual(status,200)
        self.assertNotIn(str(self.root),data.decode());self.assertNotIn(self.token,data.decode())
        self.store.logout(self.token)
        self.assertEqual(self.viewers.response(PREFIX+'/'+ticket+'/scene')[0],403)

    def docking(self,engine='vina'):
        from biomolexplorer.docking_data import write_csv
        folder=self.root/engine;folder.mkdir()
        receptor=folder/'1ABC_A.noH.pdb';receptor.write_text(atom(1,'CA','ALA','A',2,0))
        pose=folder/'pose.pdbqt';pose.write_text('REMARK VINA RESULT: -4 0 0\n'+atom(1,'C1','LIG','A',1,3))
        pdf=folder/'M1.pdf';pdf.write_bytes(b'%PDF-1.4\nfixture')
        table=folder/'docking_results.csv';write_csv(table,[dict(molecule_chembl_id='M1',canonical_smiles='CCO',
            receptor_id='1ABC_A',engine=engine,score=-4,conformer_file=pose.name,receptor_file=receptor.name,
            footprint_file=pdf.name,footprint_origin='docked_pose')])
        self.save('docking_'+engine,[table,pose,receptor,pdf]);return table,receptor,pose,pdf

    def test_docking_scene_and_pdf_are_bound_to_table_version(self):
        from biomolexplorer.docking_results import DockingResults
        table,receptor,pose,pdf=self.docking('dock6');service=DockingResults(self.store)
        version=service.page(self.token,self.pid,self.rid,self.sid,str(table))['version']
        selection=dict(table=str(table),index=0,version=version,column='conformer_file')
        key=self.viewers.issue_result(self.token,self.pid,self.rid,self.sid,'docking',selection)
        self.assertEqual(self.viewers.response(PREFIX+'/'+key+'/scene')[0],200)
        pdfkey=self.viewers.issue_document(self.token,self.pid,self.rid,self.sid,{k:v for k,v in selection.items() if k!='column'})
        response=self.viewers.response(PREFIX+'/'+pdfkey)
        self.assertEqual(response[1],'application/pdf');self.assertEqual(response[2],pdf.read_bytes())
        table.write_text(table.read_text().replace('M1','M2'))
        self.assertEqual(self.viewers.response(PREFIX+'/'+key+'/scene')[0],403)
        self.assertEqual(self.viewers.response(PREFIX+'/'+pdfkey)[0],403)

    def test_missing_receptor_does_not_silently_show_an_unrelated_structure(self):
        from biomolexplorer.docking_results import DockingResults
        table,receptor,pose,pdf=self.docking();receptor.unlink()
        service=DockingResults(self.store);version=service.page(self.token,self.pid,self.rid,self.sid,str(table))['version']
        with self.assertRaises(AccessDenied):service.scene(self.token,self.pid,self.rid,self.sid,str(table),0,version)

    def test_legacy_blank_footprint_is_not_offered_as_a_valid_graph(self):
        from biomolexplorer.docking_results import DockingResults
        table,receptor,pose,pdf=self.docking('dock6')
        plots=pdf.parent/'footprint'/'plots';plots.mkdir(parents=True)
        target=plots/pdf.name;target.write_bytes(pdf.read_bytes())
        companion=plots.parent/'M1_footprint_scored.txt';companion.write_text('')
        table.write_text(table.read_text().replace(pdf.name,'footprint/plots/'+pdf.name))
        self.save('docking_dock6',[table,receptor,pose,target,companion])
        service=DockingResults(self.store);version=service.page(self.token,self.pid,self.rid,self.sid,str(table))['version']
        with self.assertRaisesRegex(ValueError,'energias por resíduo'):
            service.footprint(self.token,self.pid,self.rid,self.sid,str(table),0,version)

    def test_legacy_pdf_is_available_even_when_manifest_omits_it(self):
        from biomolexplorer.docking_results import DockingResults
        table,receptor,pose,pdf=self.docking('dock6')
        self.save('docking_dock6',[table,receptor,pose])
        service=DockingResults(self.store)
        page=service.page(self.token,self.pid,self.rid,self.sid,str(table))
        self.assertTrue(page['rows'][0]['footprint_available'])
        self.assertEqual(service.footprint(self.token,self.pid,self.rid,self.sid,str(table),0,page['version']),str(pdf))

    def test_footprint_document_does_not_require_a_pose_or_receptor(self):
        from biomolexplorer.docking_results import DockingResults
        table,receptor,pose,pdf=self.docking('dock6')
        pose.unlink();receptor.unlink()
        self.save('consensus',[table,pdf])
        service=DockingResults(self.store)
        page=service.page(self.token,self.pid,self.rid,self.sid,str(table))
        self.assertTrue(page['rows'][0]['footprint_available'])
        key=self.viewers.issue_document(self.token,self.pid,self.rid,self.sid,
            dict(table=str(table),index=0,version=page['version']))
        self.assertEqual(self.viewers.response(PREFIX+'/'+key)[2],pdf.read_bytes())
        pdf.unlink()
        self.assertFalse(service.page(self.token,self.pid,self.rid,self.sid,str(table))['rows'][0]['footprint_available'])
