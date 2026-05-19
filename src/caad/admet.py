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
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
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
        
        # Compila os padrões SMARTS
        self.mutagenic_patterns = [Chem.MolFromSmarts(x) for x in self.MUTAGENIC_SMARTS]
        self.tumorigenic_patterns = [Chem.MolFromSmarts(x) for x in self.TUMORIGENIC_SMARTS]

    def calculate_properties(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return {
            'Mol': mol,
            'TPSA': Descriptors.TPSA(mol),
            'WLOGP': Descriptors.MolLogP(mol),
            'MW': Descriptors.MolWt(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'RB': Descriptors.NumRotatableBonds(mol)
        }

    def is_toxic(self, mol):
        """Aplica os filtros PAINS, mutagênico e tumorogênico."""
        if self.pains_catalog.GetFirstMatch(mol) is not None:
            return True
        if any(mol.HasSubstructMatch(patt) for patt in self.mutagenic_patterns):
            return True
        if any(mol.HasSubstructMatch(patt) for patt in self.tumorigenic_patterns):
            return True
        return False

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
    """Responsável exclusivamente por renderizar e salvar o gráfico BOILED-Egg."""
    
    @staticmethod
    def plot(df, output_path, output_image_file):
        path = str(Path.cwd()) + output_path + '/' + output_image_file
        fig, ax = plt.subplots(figsize=(11, 8))
        ax.set_facecolor('#d9d9d9')

        # Região HIA (Clara de ovo)
        white = Ellipse(xy=(75, 2.0), width=150, height=6.0, angle=0,
                        facecolor='white', edgecolor='black', linewidth=1.5)
        ax.add_patch(white)

        # Região BBB (Gema)
        yolk = Ellipse(xy=(42, 2.3), width=94, height=4.4, angle=0,
                       facecolor='#ffe066', edgecolor='orange', linewidth=1.5)
        ax.add_patch(yolk)

        # Plotar os pontos BBB+ e BBB-
        bbb_plus = df[df['BBB'] == 'BBB+']
        bbb_minus = df[df['BBB'] == 'BBB-']

        ax.scatter(bbb_plus['TPSA'], bbb_plus['WLOGP'], color='red', s=55,
                   edgecolors='black', linewidths=0.4, alpha=0.85, label='BBB+')
        ax.scatter(bbb_minus['TPSA'], bbb_minus['WLOGP'], color='blue', s=55,
                   edgecolors='black', linewidths=0.4, alpha=0.75, label='BBB-')

        # Configurações do gráfico
        plt.xlim(0, 200)
        plt.ylim(-2, 7)
        plt.xlabel('TPSA ($Å^2$)', fontsize=13)
        plt.ylabel('WLOGP', fontsize=13)
        plt.title('BOILED-Egg Plot', fontsize=15)
        plt.legend(loc='upper right')
        plt.tight_layout()

        plt.savefig(path, dpi=300)
        plt.close() # Libera memória do pyplot