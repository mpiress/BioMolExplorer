"""Bounded, cancellable structural queries independent of the dashboard controls."""
from biomolexplorer.molecular_parsing import mol_from_smiles
import base64
import io


def depiction(molecule, query=None, size=(400,280)):
    from rdkit.Chem import Draw
    atoms=list(molecule.GetSubstructMatch(query)) if query is not None else []
    bonds=[molecule.GetBondBetweenAtoms(atoms[b.GetBeginAtomIdx()],atoms[b.GetEndAtomIdx()]).GetIdx() for b in query.GetBonds()] if atoms else []
    out=io.BytesIO()
    Draw.MolToImage(molecule,size=size,highlightAtoms=atoms,highlightBonds=bonds).save(out,format='PNG')
    return out.getvalue()


def individual(smiles, smarts=''):
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:return dict(image=None,scaffold=None,atoms=0,smiles='',notice='RDKit indisponível neste ambiente.')
    molecule=mol_from_smiles(smiles) if isinstance(smiles,str) else None
    if molecule is None:return dict(image=None,scaffold=None,atoms=0,smiles='',notice='SMILES ausente ou inválido.')
    scaffold=MurckoScaffold.GetScaffoldForMol(molecule)
    query=Chem.MolFromSmarts(smarts) if smarts else None
    match=molecule.GetSubstructMatch(query) if query is not None else ()
    return dict(image=depiction(molecule,query if match else None),
        scaffold=depiction(scaffold) if scaffold.GetNumAtoms() else None,
        atoms=scaffold.GetNumAtoms(),smiles=Chem.MolToSmiles(scaffold),
        notice='' if scaffold.GetNumAtoms() else 'Molécula acíclica: scaffold Bemis–Murcko vazio.',
        highlighted_atoms=len(match))


def common(rows, timeout, limit, cancelled):
    """Threshold 1.0: every molecule in the explicit set must match; never sample."""
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    result=dict(ids=[n for n,_ in rows],compounds=len(rows),timeout=timeout,limit=limit,
        ring_matches_ring_only=True,complete_rings_only=True,threshold=1.0,atom_compare='elements',bond_compare='order',match_valences=False,
        smarts='',image='',atoms=0,bonds=0,status='unavailable')
    if not 1<=timeout<=30 or limit<2:
        return dict(result,notice='Use 1–30 segundos e limite de pelo menos 2 compostos.')
    if len(rows)<2 or len(rows)>limit:
        return dict(result,notice='Selecione de 2 até o limite de compostos. Nenhum subconjunto é amostrado automaticamente.')
    molecules=[]
    for _,smiles in rows:
        if cancelled.is_set():return dict(result,status='cancelled',notice='Cálculo interrompido.')
        molecule=mol_from_smiles(smiles) if isinstance(smiles,str) else None
        if molecule is None:return dict(result,notice='Há SMILES ausentes ou inválidos no conjunto; nenhum composto foi excluído silenciosamente.')
        if molecule.GetNumAtoms()>256:return dict(result,notice='Estrutura acima de 256 átomos: restrinja o conjunto para esta análise interativa.')
        molecules.append(molecule)
    class Progress(rdFMCS.MCSProgress):
        def __call__(self,progress,parameters):return not cancelled.is_set()
    parameters=rdFMCS.MCSParameters();parameters.Timeout=timeout;parameters.Threshold=1.0
    parameters.AtomCompareParameters.RingMatchesRingOnly=True
    parameters.BondCompareParameters.RingMatchesRingOnly=True
    parameters.BondCompareParameters.CompleteRingsOnly=True
    parameters.ProgressCallback=Progress()
    found=rdFMCS.FindMCS(molecules,parameters)
    if cancelled.is_set():return dict(result,status='cancelled',notice='Cálculo interrompido.')
    query=Chem.MolFromSmarts(found.smartsString) if found.smartsString else None
    if query is None or found.numAtoms<2:
        return dict(result,status='partial' if found.canceled else 'no_common_fragment',notice='Nenhuma subestrutura comum informativa (≥ 2 átomos). Explore um conjunto menor.')
    matches=sum(bool(mol.GetSubstructMatch(query)) for mol in molecules)
    if matches!=len(molecules):return dict(result,notice='Consulta sem cobertura completa; resultado não apresentado como MCS do conjunto.')
    from caad.graph_results import fragment_molecule
    fragment=fragment_molecule(found.smartsString,rows[0][1])
    return dict(result,status='partial' if found.canceled else 'complete',smarts=found.smartsString,
        atoms=found.numAtoms,bonds=found.numBonds,smiles=Chem.MolToSmiles(fragment),image=base64.b64encode(depiction(fragment)).decode('ascii'),
        coverage=matches/len(molecules),reference_id=rows[0][0],notice='Tempo limite: melhor fragmento encontrado; maximalidade não confirmada.' if found.canceled else 'Busca concluída; cobertura de 100% do conjunto explícito.')
