from biomolexplorer.paths import directory, resolve_path, worker_count
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

   module_description=(
    "Wrapper module for ADMET-based molecular evaluation "
    "functionalities available in the src/caad diretory"
),

    module_version="1.0.0"
)

#----------------------------------------------------------------------------------------------
from pathlib import Path
import pandas as pd
import os
import json
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from caad.admet import MoleculeEvaluator, BoiledEggPlotter
from kernel.utilities import fileHandling
from biomolexplorer.visualizations import egg_view, write_view, SUFFIX
#----------------------------------------------------------------------------------------------

class ADMETWrapper:

    def __init__(self, base_output_path, base_input_path=None, input_file=None, verbose=0,
                 admet_profile='legacy', admet_mode=None, prediction_backend='rdkit',
                 bbb_target=None, hia_target=None, rule_set=None, max_rule_violations=None,
                 alert_policy=None, alert_catalogs=None, property_limits=None,dataset_origin=None):
        from biomolexplorer.admet_config import resolve_config
        self.config=resolve_config(dict(admet_profile=admet_profile,admet_mode=admet_mode,
            prediction_backend=prediction_backend,bbb_target=bbb_target,hia_target=hia_target,
            rule_set=rule_set,max_rule_violations=max_rule_violations,alert_policy=alert_policy,
            alert_catalogs=alert_catalogs,property_limits=property_limits))
        self.base_output_path = base_output_path
        self.set_output_path()

        self.base_input_path = base_input_path
        self.dataset_origin = dataset_origin
        self.csv_files = self.set_input_path(input_file)

        self.evaluator = MoleculeEvaluator()
        from biomolexplorer.calculation_cache import CalculationCache
        from biomolexplorer.stage_cache import file_digest
        from rdkit import rdBase
        import inspect
        self.structure_cache=CalculationCache('admet-descriptors-v1',self.base_output_path)
        self.structure_recipe=dict(rdkit=rdBase.rdkitVersion,
            evaluator=file_digest(inspect.getfile(MoleculeEvaluator)),wrapper=file_digest(__file__),
            alerts=self.config['alert_catalogs'],legacy=self.config['admet_profile']=='legacy')
        self.structure_results={}
        self.df_results = None
        self.excluded_count = 0

        self.verbose = 0


    def set_input_path(self, input_file) -> list:

        base_input_path = Path(directory(self.base_input_path))

        if not base_input_path.is_dir():
            raise FileNotFoundError(f"O caminho especificado não existe ou não é um diretório: {self.base_input_path}")

        csv_files = [str(arquivo.name) for arquivo in base_input_path.glob('*.csv')] if input_file is None else [input_file]
        return csv_files


    def set_output_path(self) -> None:
        base_output_path = Path(directory(self.base_output_path))
        if not base_output_path.exists() or not self.base_output_path is None:
           os.makedirs(base_output_path, exist_ok=True)


    def _load_data(self, filename):
        f1 = fileHandling(input_path=self.base_input_path, output_path=self.base_output_path)
        df = f1.csv_to_dataframe(filename)

        missing = {'molecule_chembl_id','canonical_smiles'} - set(df.columns)
        if missing:
            raise ValueError('O CSV de compostos precisa das colunas: ' + ', '.join(sorted(missing)))

        from biomolexplorer.molecule_quality import clean_dataframe
        return clean_dataframe(df,source=str(Path(self.base_input_path)/(filename+'.csv')))


    def run_pipeline(self):
        from biomolexplorer.admet_config import selection_reasons,DESCRIPTORS
        from rdkit import rdBase
        self.total_evaluated=0;self.total_hypothesis_selected=0
        for filename in self.csv_files:
            filename=Path(filename).stem;df_input=self._load_data(filename)
            from biomolexplorer.dataset_identity import output_stem,admet_origin
            origin=self.dataset_origin or admet_origin([Path(self.base_input_path)/(filename+'.csv')])
            filename=output_stem(filename,origin)
            rows=[];self.excluded_count=0
            provider=dict(provider='RDKit',version=rdBase.rdkitVersion,domain_assessment='not_assessed')
            predictions={}
            if self.config['prediction_backend']=='admet_ai' and len(df_input):
                from biomolexplorer.admet_predictors import predict_admet_ai
                predictions,provider=predict_admet_ai(df_input['canonical_smiles'].tolist())
            for _,row in df_input.iterrows():
                smiles=row['canonical_smiles'];properties=self.evaluate_structure(smiles)
                if properties is None:continue
                legacy_toxic=properties.pop('_legacy_toxic')
                result=row.to_dict()
                result.update(properties)
                result.update(PGP=properties['PGP'] if self.config['admet_profile']=='legacy' else 'not_evaluated',
                    BBB_method='BioMolExplorer legacy ellipse + penalties v1',HIA_method='BioMolExplorer local ellipse v1',
                    ADMET_profile=self.config['admet_profile'],Prediction_backend=self.config['prediction_backend'])
                result.update(predictions.get(smiles,{}))
                reasons=selection_reasons(result,self.config)
                if self.config['admet_profile']=='legacy':
                    result['PGP_method']='Legacy descriptor heuristic; not a validated substrate model'
                    if legacy_toxic:reasons.append('Legacy structural exclusion (not measured toxicity)')
                result.update(ADMET_selected=not reasons,ADMET_reasons='; '.join(reasons))
                rows.append(result)
            columns=list(dict.fromkeys([*df_input.columns,*DESCRIPTORS,'BBB','HIA','PGP','BBB_method','HIA_method',
                'StructuralAlert_count','StructuralAlerts','ADMET_warnings','ADMET_profile','Prediction_backend','ADMET_selected','ADMET_reasons',
                *['Alert_'+name for name in self.config['alert_catalogs']],
                *[key for row in rows for key in row]]))
            evaluated=pd.DataFrame(rows,columns=columns)
            selected=evaluated[evaluated['ADMET_selected']==True].copy()
            self.excluded_count=len(evaluated)-len(selected)
            self.total_evaluated+=len(evaluated);self.total_hypothesis_selected+=len(selected)
            results=selected if self.config['admet_mode']=='filter' else evaluated
            self.df_results=results;self.evaluated=evaluated
            # Preserve every valid molecule and each failed criterion even when filtering.
            evaluated.to_csv(Path(self.base_output_path)/(filename+'_evaluated.csv'),index=False)
            self.export_results(results,filename);self.generate_plot(results,filename)
            from biomolexplorer.admet_outputs import export_standard,SCHEMA_VERSION
            outputs=export_standard(evaluated,self.config,self.base_output_path,filename)
            metadata=dict(schema_version=SCHEMA_VERSION,outputs=outputs,dataset_origin=origin,output_stem=filename,
                configuration=self.config,provider=provider,rdkit_version=rdBase.rdkitVersion,
                evaluated=len(evaluated),selected=len(selected),output=len(results),
                descriptors=DESCRIPTORS,methods=dict(
                    BBB='Local legacy ellipse with penalties; not validated BOILED-Egg',
                    HIA='Local ellipse approximation; passive permeability only',
                    ESOL_logS='RDKit ESOL-style refit by Pat Walters, log10 mol/L',
                    logKp='Potts–Guy-style regression with RDKit logP substitution, log10 cm/s',
                    alerts='RDKit PAINS/BRENK/NIH/ZINC/ChEMBL matches; not proof of toxicity'),
                unavailable_without_trained_model=['CYP inhibition/substrate','hERG','AMES','DILI','clearance','half-life','PPBR','VDss'],
                limitations=['No pH-dependent logD/pKa or unbound brain exposure model',
                    'No validated applicability-domain or prediction-interval assessment',
                    'Rules and thresholds are research hypotheses, not efficacy/safety certification'])
            (Path(self.base_output_path)/(filename+'_admet_report.json')).write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n')
        return self.df_results

    def evaluate_structure(self,smiles):
        """Reuse molecular descriptors across subsets, merges and pipeline runs.

        Input annotations and selection thresholds are deliberately applied
        outside this cache; they must always reflect the current configuration.
        """
        from biomolexplorer.calculation_cache import identity
        key=identity(dict(self.structure_recipe,smiles=smiles))
        cached=self.structure_results.get(key) or self.structure_cache.metadata(key)
        if cached and 'properties' in cached:return dict(cached['properties'])
        props=self.evaluator.calculate_properties(smiles)
        if props is None:return None
        result={k:v for k,v in props.items() if k!='Mol'}
        result.update(self.evaluator.screening_properties(props,self.config['alert_catalogs']))
        result.update(BBB=self.evaluator.classify_bbb(props['TPSA'],props['WLOGP'],props['MW'],props['HBD'],props['RB']),
            HIA=self.evaluator.classify_hia(props['TPSA'],props['WLOGP']),
            PGP=self.evaluator.predict_pgp(props['MW'],props['TPSA'],props['WLOGP'],props['HBD']) if self.structure_recipe['legacy'] else 'not_evaluated',
            _legacy_toxic=bool(self.evaluator.is_toxic(props['Mol'])) if self.structure_recipe['legacy'] else False)
        cached=dict(properties=result)
        self.structure_cache.remember(key,{},cached)
        self.structure_results[key]=cached
        return dict(result)

    def export_results(self, evaluated, filename):
        if evaluated is None:raise ValueError('Execute a avaliação ADMET antes de exportar.')
        self.output_files={}
        def write(role,frame):
            path=Path(self.base_output_path)/(filename+'_'+role+'.csv')
            frame.to_csv(path,index=False);self.output_files[role]=str(path)
        for field,label in [('BBB','BBB+'),('BBB','BBB-'),('HIA','HIA+'),('HIA','HIA-')]:
            write(label,evaluated[evaluated[field]==label])
        rule=self.config['rule_set']
        if rule!='none':
            names={'lipinski':['Lipinski'],'veber':['Veber'],'oral':['Lipinski','Veber'],'fragment':['Fragment']}[rule]
            passed=pd.Series(True,index=evaluated.index,dtype=bool)
            for name in names:
                allowed=0 if rule=='oral' and name=='Veber' else self.config['max_rule_violations']
                passed &= evaluated[name+'_violations']<=allowed
            rule_label='lipinski-veber' if rule=='oral' else rule
            table=evaluated[passed].copy();table['Rule_set']=rule_label
            write('rules_'+rule_label,table)
        if self.config['alert_catalogs']:
            table=evaluated.copy()
            for catalog in self.config['alert_catalogs']:
                column='Alert_'+catalog
                table[column]=table[column].map(lambda value:'true' if value is True or str(value).lower()=='true' else 'false')
            write('alerts',table)


    def generate_plot(self, results, filename):
        if results is None:
            raise ValueError("Nenhum dado processado para gerar o gráfico.")
        BoiledEggPlotter.plot(results, self.base_output_path, filename+'_egg.png')
        write_view(Path(self.base_output_path) / (filename + '_egg' + SUFFIX),
                   egg_view(results, filename + ' · BOILED-Egg'))


    def print_summary(self, results):
        if results is None:
            return
        total_valid = len(results)
        print('\n========== SUMMARY ==========')
        print(f'Total molecules input: {total_valid + self.excluded_count}')
        print(f'Compounds outside configured research criteria: {self.excluded_count}')
        print(f'Molecules processed: {total_valid}')
        print(f'BBB+: {len(results[results["BBB"]=="BBB+"])}')
        print(f'BBB-: {len(results[results["BBB"]=="BBB-"])}')
        print(f'HIA+: {len(results[results["HIA"]=="HIA+"])}')
        print('=============================\n')
