from biomolexplorer.molecular_parsing import mol_from_smiles
from biomolexplorer.paths import directory
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
        "ADMET molecule evaluator"
    ),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
from pathlib import Path
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.loggers import LoggerManager
#----------------------------------------------------------------------------------------------


class MoleculeEvaluator:
    """Responsável por calcular propriedades e aplicar filtros em moléculas (SMILES)."""

    MUTAGENIC_SMARTS = ['[N+](=O)[O-]', 'N=N', '[CX3](=O)[Cl]', '[SH]', '[C,c]Br', '[C,c]I']
    TUMORIGENIC_SMARTS = ['[N+](=O)[O-]', 'C=C=O', '[Cl][C]=O', '[C,c]Cl', '[C,c]Br']

    def __init__(self):
        # Inicializa o catálogo PAINS
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        self.pains_catalog = FilterCatalog(params)
        self.alert_catalogs = {'PAINS':self.pains_catalog}
        for name in ('BRENK','NIH'):
            params=FilterCatalogParams();params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs,name))
            self.alert_catalogs[name]=FilterCatalog(params)

        # Compila os padrões SMARTS
        self.mutagenic_patterns = [Chem.MolFromSmarts(x) for x in self.MUTAGENIC_SMARTS]
        self.tumorigenic_patterns = [Chem.MolFromSmarts(x) for x in self.TUMORIGENIC_SMARTS]

        self.logger = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/admet.log')


    def calculate_properties(self, smiles):

        try:
            mol = mol_from_smiles(smiles)
            if mol is None:
                return None
            return {
                'Mol': mol,
                'TPSA': Descriptors.TPSA(mol),
                'WLOGP': Descriptors.MolLogP(mol),
                'MW': Descriptors.MolWt(mol),
                'HBD': Descriptors.NumHDonors(mol),
                'HBA': Descriptors.NumHAcceptors(mol),
                'RB': Descriptors.NumRotatableBonds(mol),
                'QED': QED.qed(mol),
                'FractionCSP3': rdMolDescriptors.CalcFractionCSP3(mol),
                'HeavyAtoms': mol.GetNumHeavyAtoms(),
                'Rings': rdMolDescriptors.CalcNumRings(mol),
                'AromaticRings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'FormalCharge': Chem.GetFormalCharge(mol),
                'Fragments': len(Chem.GetMolFrags(mol)),
                'AromaticProportion': sum(a.GetIsAromatic() for a in mol.GetAtoms())/max(1,mol.GetNumHeavyAtoms())
            }

        except Exception as e:
            self.logger.error(f'Error during to perform the calculate_properties function', exc_info=True)
            raise


    def screening_properties(self, props, catalogs):
        """Descriptors and interpretable rules; alerts are not toxicity labels."""
        mol=props['Mol'];mw,logp,tpsa=props['MW'],props['WLOGP'],props['TPSA']
        hbd,hba,rb=props['HBD'],props['HBA'],props['RB']
        # ESOL-style RDKit refit documented by Pat Walters (see docs/admet_research.md).
        values=dict(
            ESOL_logS=0.26121066137801696-0.0066138847738667125*mw-0.7416739523408995*logp+
                0.003451545565957996*rb-0.42624840441316975*props['AromaticProportion'],
            logKp=-6.3+0.71*logp-0.0061*mw,  # Potts–Guy-style, RDKit logP substitution (cm/s).
            Lipinski_violations=sum((mw>500,logp>5,hbd>5,hba>10)),
            Veber_violations=sum((rb>10,tpsa>140)),
            Fragment_violations=sum((mw>300,logp>3,hbd>3,hba>3,rb>3,tpsa>60)))
        alerts=[name+': '+entry.GetDescription() for name in catalogs for entry in self.alert_catalogs[name].GetMatches(mol)]
        values.update(StructuralAlert_count=len(alerts),StructuralAlerts='; '.join(alerts))
        warnings=[]
        if props['Fragments']>1:warnings.append('Disconnected structure: salt/mixture; descriptors describe the submitted structure')
        if any(atom.GetFormalCharge() for atom in mol.GetAtoms()):warnings.append('Charged structure: pH-dependent permeability is not modeled')
        if mw>700 or tpsa>200:warnings.append('Large/polar structure: passive-permeability heuristics require particular caution')
        values['ADMET_warnings']='; '.join(warnings)
        return values

    def is_toxic(self, mol):

        try:
            if self.pains_catalog.GetFirstMatch(mol) is not None:
                return True
            if any(mol.HasSubstructMatch(patt) for patt in self.mutagenic_patterns):
                return True
            if any(mol.HasSubstructMatch(patt) for patt in self.tumorigenic_patterns):
                return True
            return False

        except Exception as e:
            self.logger.error(f'Error during to perform the is_toxic function', exc_info=True)
            raise


    def predict_pgp(self, mw, tpsa, logp, hbd):
        score = sum([mw > 400, tpsa > 75, hbd >= 2, logp > 3.5])
        return "PGP+" if score >= 2 else "PGP-"

    def classify_bbb(self, tpsa, logp, mw, hbd, rb):
        # Elipse interna do BOILED-Egg para BBB
        ellipse_value = (((tpsa - 42) ** 2) / (47 ** 2)) + (((logp - 2.3) ** 2) / (2.2 ** 2))
        inside_ellipse = ellipse_value <= 1

        # Penalidades heurísticas
        penalty = sum([tpsa > 90, mw > 500, rb > 10, hbd > 2, logp < 0, logp > 6])
        return "BBB+" if (inside_ellipse and penalty <= 1) else "BBB-"

    def classify_hia(self, tpsa, logp):
        value = (((tpsa - 75) ** 2) / (75 ** 2)) + (((logp - 2.0) ** 2) / (3.0 ** 2))
        return "HIA+" if value <= 1 else "HIA-"


class BoiledEggPlotter:
    @staticmethod
    def plot(df, output_path, output_image_file):
        from biomolexplorer.egg_plot import egg_png
        from biomolexplorer.visualizations import egg_view
        Path(directory(output_path),output_image_file).write_bytes(egg_png(egg_view(df,'EGG · BioMolExplorer')))
