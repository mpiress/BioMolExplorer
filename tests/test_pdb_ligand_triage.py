"""Geometry-driven PDB shortlists retain alternatives and scientific evidence."""
import asyncio
import csv
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import flet as ft
from biomolexplorer.catalog import new_stage
from biomolexplorer.operations import validate_operation
from biomolexplorer.pdb_ligands import analyze_ligands
from biomolexplorer.ui.guided import GuidedForm
from crawlers.complex import PDBComplex
import test_retrieval_usability as usability


def atom(serial,name,residue,chain,number,xyz,element='C',het=False,icode='',occupancy=1):
    x,y,z=xyz
    return (f"{'HETATM' if het else 'ATOM':6s}{serial:5d} {name:^4s} {residue:>3s} {chain}{number:4d}{icode:1s}   "
            f'{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{20.:6.2f}          {element:>2s}  \n')


def pocket(chain='A',offset=0,start=1,serial=1):
    result=''
    for i,(name,xyz) in enumerate(zip(('ALA','GLY','SER','VAL'),((-2,0,0),(2,0,0),(0,-2,0),(0,2,0)))):
        result+=atom(serial+i,'CA',name,chain,start+i,(xyz[0]+offset,xyz[1],xyz[2]))
    return result


def ligand(code='LIG',chain='A',number=101,offset=0,serial=100,**kwargs):
    return ''.join(atom(serial+i,f'C{i+1}',code,chain,number,(offset+i*.15,0,.5),het=True,**kwargs) for i in range(8))


def site_record(site,residues):
    line=list(' '*80);line[:6]='SITE  ';line[7:10]=f'{1:3d}';line[11:14]=f'{site:3s}';line[15:17]=f'{len(residues):2d}'
    for offset,(name,chain,number) in zip((18,29,40,51),residues):
        line[offset:offset+3]=f'{name:>3s}';line[offset+4]=chain;line[offset+5:offset+9]=f'{number:4d}'
    return ''.join(line)+'\n'


def link_record(left,right):
    line=list(' '*80);line[:6]='LINK  '
    for offset,(name,chain,number) in zip((17,47),(left,right)):
        line[offset:offset+3]=f'{name:>3s}';line[offset+4]=chain;line[offset+5:offset+9]=f'{number:4d}'
    return ''.join(line)+'\n'


def rows(text):return analyze_ligands(text,'1ABC')['candidates']


class LigandTriageTests(unittest.TestCase):
    def test_prioritizes_bound_organic_ligand_and_explains_exclusions(self):
        text=pocket()+ligand()+ligand('GOL',number=102,serial=110)
        text+=atom(120,'ZN','ZN','A',103,(0,0,1),'ZN',het=True)
        text+=ligand('FAR',number=104,offset=50,serial=130)
        text+=atom(140,'O','HOH','A',105,(0,0,1),'O',het=True)
        result=rows(text);by_name={r['LIGAND']:r for r in result}
        self.assertEqual([r['LIGAND'] for r in result if r['RECOMMENDED']],['LIG'])
        self.assertEqual(by_name['GOL']['REASON'],'additive')
        self.assertEqual(by_name['ZN']['REASON'],'nonorganic')
        self.assertEqual(by_name['FAR']['REASON'],'contacts')
        self.assertNotIn('HOH',by_name)
        self.assertEqual(by_name['LIG']['CONTACT_RESIDUES'],4)
        self.assertEqual(by_name['LIG']['CONTACT_FRACTION'],1)
        self.assertEqual(by_name['LIG']['ACTIVE_SITE_CONTACTS'],0)

    def test_deduplicates_equivalent_copies_but_retains_distinct_pockets(self):
        text=pocket()+ligand()+pocket('B',20,serial=20)+ligand(chain='B',offset=20,serial=120)
        text+=pocket('A',50,start=11,serial=30)+ligand(number=102,offset=50,serial=140)
        text+=pocket('B',70,start=11,serial=50)
        result=rows(text)
        selected=[r for r in result if r['RECOMMENDED']]
        self.assertEqual(len(selected),2)
        duplicate=next(r for r in result if r['REASON']=='equivalent_copy')
        self.assertEqual(duplicate['CHAIN'],'B')
        self.assertEqual(duplicate['REPRESENTATIVE'],'LIG:A:101')

    def test_site_annotation_prioritizes_catalytic_pocket_without_inventing_activity(self):
        annotation='REMARK 800 SITE_IDENTIFIER: AC1\nREMARK 800 SITE_DESCRIPTION: CATALYTIC SITE\n'
        annotation+=site_record('AC1',[('ALA','A',1),('GLY','A',2),('SER','A',3)])
        text=annotation+pocket()+ligand()+pocket(offset=25,start=11,serial=20)+ligand('DRG',number=102,offset=25,serial=120)
        result=rows(text);selected=[r for r in result if r['RECOMMENDED']]
        self.assertEqual([r['LIGAND'] for r in selected],['LIG'])
        self.assertEqual(selected[0]['REASON'],'active_site')
        self.assertEqual(selected[0]['SITE_IDS'],'AC1')
        self.assertEqual(next(r for r in result if r['LIGAND']=='DRG')['REASON'],'other_site')
        generic=rows(text.replace('CATALYTIC SITE','BINDING SITE FOR RESIDUE LIG A 101'))
        self.assertEqual(len([r for r in generic if r['RECOMMENDED']]),2)
        self.assertTrue(all(r['ACTIVE_SITE_CONTACTS']==0 for r in generic))

    def test_covalent_link_and_conect_do_not_enter_conventional_redocking_shortlist(self):
        text=pocket()+ligand()
        linked=rows(link_record(('ALA','A',1),('LIG','A',101))+text)[0]
        connected=rows(text+'CONECT    1  100\n')[0]
        for row in (linked,connected):
            self.assertTrue(row['PROTEIN_LINKED']);self.assertFalse(row['RECOMMENDED'])
            self.assertEqual(row['REASON'],'protein_linked')

    def test_cofactor_remains_an_option_and_is_a_fallback_when_alone(self):
        text=pocket()+ligand('FAD')
        self.assertTrue(rows(text)[0]['RECOMMENDED'])
        result=rows(text+ligand('DRG',number=102,serial=120))
        self.assertEqual([r['LIGAND'] for r in result if r['RECOMMENDED']],['DRG'])
        self.assertEqual(next(r for r in result if r['LIGAND']=='FAD')['REASON'],'cofactor')
        # ATP can be a substrate/inhibitor; do not blacklist nucleotides.
        self.assertTrue(rows(pocket()+ligand('ATP'))[0]['RECOMMENDED'])

    def test_linked_glycan_branch_is_not_mistaken_for_a_free_ligand(self):
        text=pocket()+ligand('NAG')+ligand('MAN',number=102,serial=120)
        text+=link_record(('ALA','A',1),('NAG','A',101))
        text+=link_record(('NAG','A',101),('MAN','A',102))
        result=rows(text)
        self.assertTrue(all(r['REASON']=='protein_linked' for r in result))

    def test_contact_fraction_and_distinct_residue_thresholds(self):
        one_residue=atom(1,'CA','ALA','A',1,(0,0,0))+ligand()
        self.assertEqual(rows(one_residue)[0]['REASON'],'contacts')
        text=pocket()+ligand()
        # Most of this residue's observed atoms are away from the pocket.
        text+=''.join(atom(200+i,f'N{i}','LIG','A',101,(50+i,0,0),'N',het=True) for i in range(12))
        row=rows(text)[0]
        self.assertEqual(row['CONTACT_RESIDUES'],4)
        self.assertLess(row['CONTACT_FRACTION'],.5)
        self.assertFalse(row['RECOMMENDED'])

    def test_first_model_zero_occupancy_hydrogens_and_insertion_codes(self):
        first=pocket()+ligand(offset=50)
        second=pocket()+ligand()
        result=rows('MODEL        1\n'+first+'ENDMDL\nMODEL        2\n'+second+'ENDMDL\n')
        self.assertFalse(result[0]['RECOMMENDED'])
        result=rows(pocket()+ligand(occupancy=0))
        self.assertEqual(result[0]['HEAVY_ATOMS'],0)
        result=rows(pocket()+ligand(element='H'))
        self.assertFalse(result[0]['RECOMMENDED'])
        result=rows(pocket()+ligand(icode='A')+ligand(icode='B',serial=120))
        self.assertEqual(len(result),2)
        self.assertTrue(all(r['REASON']=='unsupported_identity' for r in result))

    def test_different_protein_chain_is_reported_without_silently_relabeling_ligand(self):
        result=rows(pocket()+ligand(chain='B'))[0]
        self.assertEqual(result['CHAIN'],'B')
        self.assertEqual(result['RECEPTOR_CHAINS'],'A')
        self.assertEqual(result['REASON'],'different_chain')
        self.assertFalse(result['RECOMMENDED'])

    def test_empty_shortlist_is_explicit_and_keeps_audit_rows(self):
        result=rows(pocket()+ligand('EDO'))
        self.assertFalse(any(r['RECOMMENDED'] for r in result))
        self.assertEqual(len(result),1)

    def test_deposited_component_name_is_retained_for_manual_review(self):
        text='HETNAM     LIG EXAMPLE ORGANIC\nHETNAM   2 LIG COMPOUND\n'+pocket()+ligand()
        self.assertEqual(rows(text)[0]['LIGAND_NAME'],'EXAMPLE ORGANIC COMPOUND')

    def test_one_kdtree_is_built_for_each_structure(self):
        from Bio.PDB import NeighborSearch
        with patch('biomolexplorer.pdb_ligands.NeighborSearch',wraps=NeighborSearch) as tree:
            rows(pocket()+ligand()+ligand('DRG',number=102,serial=120))
        tree.assert_called_once()

    def test_retrieval_writes_shortlist_full_evidence_and_empty_structure_notice(self):
        with tempfile.TemporaryDirectory() as folder:
            crawler=PDBComplex(folder);crawler.get_pdb_ids_with_filters=Mock(return_value=['1ABC','2ABC'])
            text=pocket()+ligand()+ligand('GOL',number=102,serial=120)
            session=Mock()
            session.get.side_effect=lambda url,**kwargs:SimpleNamespace(
                text=text if '1ABC' in url else pocket()+ligand('GOL'),raise_for_status=lambda:None)
            with patch('requests.Session',return_value=session):report=crawler.get_pdb_files({'pdb_ids':['1ABC','2ABC']})
            root=Path(folder)
            selected=list(csv.DictReader(io.StringIO((root/'pdb_codes.csv').read_text())))
            all_rows=list(csv.DictReader(io.StringIO((root/'ligand_candidates.csv').read_text())))
            self.assertEqual([(r['PDB_CODE'],r['LIGAND']) for r in selected],[('1ABC','LIG')])
            self.assertEqual(len(all_rows),3)
            self.assertEqual(report['ligand_triage']['structures_without_recommendation'],['2ABC'])
            self.assertEqual(report['ligand_triage']['selected'],1)
            self.assertEqual(json.loads((root/'retrieval_report.json').read_text())['ligand_triage']['contact_distance_angstrom'],4)
            self.assertEqual(session.get.call_count,2)  # No per-ligand provider requests.

    def test_all_mode_and_parameter_validation(self):
        validate_operation('retrieve_structures',{'pdb_ids':['1ABC'],'ligand_selection':'all'})
        with self.assertRaises(ValueError):validate_operation('retrieve_structures',{'pdb_ids':['1ABC'],'ligand_selection':'bad'})
        stage=new_stage('retrieve_structures')
        self.assertEqual(stage['parameters']['ligand_selection'],'recommended')
        form=GuidedForm(SimpleNamespace(page=SimpleNamespace(update=lambda:None)),stage,[],True)
        self.assertEqual(form.field_controls['ligand_selection'].value,'recommended')
        self.assertEqual({o.key for o in form.field_controls['ligand_selection'].options},{'all','recommended'})


class TriageCurationTests(unittest.TestCase):
    """Exercise the authorized curation and pipeline handoff with realistic atoms."""
    def setUp(self):
        usability.PDBCurationTests.setUp(self)
        self.path.write_text(pocket()+ligand()+ligand('GOL',number=102,serial=120))
        self.csv.write_text('PDB_CODE,LIGAND,RESNUM,CHAIN,RESOLUTION\n1ABC,LIG,101,A,2.0\n')

    tearDown=usability.PDBCurationTests.tearDown
    context=usability.PDBCurationTests.context
    save=usability.PDBCurationTests.save
    def test_recommendation_can_be_replaced_by_an_excluded_but_present_ligand(self):
        context=self.context()
        self.assertEqual([r['LIGAND'] for r in context['recommended_ligands']],['LIG'])
        self.assertEqual({r['LIGAND'] for r in context['available_ligands']},{'LIG','GOL'})
        updated=self.save([{'LIGAND':'GOL','RESNUM':102,'CHAIN':'A'}])
        self.assertEqual([r['LIGAND'] for r in updated['ligands']],['GOL'])
        from biomolexplorer.redocking_config import metadata_records
        self.assertEqual(metadata_records([self.csv])[0][:4],['1ABC','GOL',102,'A'])

    def test_dialog_filters_options_and_restores_recommended_without_erasing_other_choices(self):
        from biomolexplorer.ui.pdb_results import PDBActions
        dialogs=[]
        async def call(function,*args):return function(*args)
        async def guard(action):await action()
        ui=SimpleNamespace(call=call,guard=guard,notify=Mock(),editing_stage=False,
            page=SimpleNamespace(update=lambda:None,show_dialog=dialogs.append))
        results=SimpleNamespace(ui=ui,token=self.token,project_id=self.pid,run_id=self.rid,stage=self.stage,service=self.service,valid=lambda:True)
        asyncio.run(PDBActions(results,True).open_ligands({'path':str(self.path)}))
        controls=dialogs[-1].content.content.controls
        listing,choice=controls[1:3];toggle=controls[4]
        self.assertEqual(len(choice.options),1)
        self.assertIn('LIG',choice.options[0].text)
        toggle.value=True;toggle.on_change(None)
        self.assertEqual(len(choice.options),2)
        self.assertEqual(len(listing.controls),1)
        choice.value='0';controls[3].on_click(None)  # Alphabetical first: GOL.
        self.assertEqual(len(listing.controls),2)
        controls[-1].on_click(None)
        self.assertEqual(len(listing.controls),1)
        self.assertEqual(listing.controls[0].controls[0].value,'LIG')

    def test_curation_does_not_discard_an_insertion_code_and_select_a_different_residue(self):
        self.path.write_text(pocket()+ligand(icode='A'))
        with self.assertRaisesRegex(ValueError,'código de inserção'):
            self.save([{'LIGAND':'LIG','RESNUM':101,'CHAIN':'A'}])


if __name__=='__main__':unittest.main()
