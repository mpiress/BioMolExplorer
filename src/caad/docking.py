from biomolexplorer.molecular_surface import generate_surface
from biomolexplorer.processes import run_command
from biomolexplorer.storage import sync_directory
from biomolexplorer.paths import directory, resolve_path, worker_count
#----------------------------------------------------------------------------------------------
#Configure PYTHONPATH to perform execution using the project classes
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="ADMET analysis",

    module_description=(
    "Module for managing and integrating molecular "
    "analysis by consensus docking strategy and redocking by AutoDock Vina"
),

    module_version="1.0.0"
)
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
import warnings
# Desabilita todos os avisos

import matplotlib
matplotlib.use('agg')
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
import subprocess
import os
import math
import time
import matplotlib.pyplot as plt
import numpy as np

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from typing import Optional, List, Tuple, Literal
from pathlib import Path
from pandas import DataFrame
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from pymol import cmd
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
from kernel.utilities import fileHandling, MolConverter, MolExplorer
from kernel.loggers import LoggerManager
from biomolexplorer.diagnostics import log_context, diagnose_exception, event
from kernel.descriptors import Descriptors
#----------------------------------------------------------------------------------------------


class Docking():

    def __init__(self, ligand_input_path:Optional[str]=None, receptor_input_path:Optional[str]=None,
                 complex_input_path:Optional[str]=None, output_path:Optional[str]=None,
                 mol_filename:Optional[str]='molecules') -> None:

        self.path           = str(Path.cwd())
        self.ligandpath     = ligand_input_path if ligand_input_path != None else None
        self.receptorpath   = receptor_input_path if receptor_input_path != None else None
        self.complexpath    = complex_input_path if complex_input_path != None else None
        self.nprocess       = worker_count()
        self.centers        = None
        self.logpath        = self.path + '/logs/'
        self.logger         = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/docking.log')
        self.mol_filename   = mol_filename


        self.set_outputpath(output_path) if output_path != None else None




    def set_ligandpath(self, path) -> None:
       if not os.path.exists(directory(path)):
           self.logger.error('Missing ligand directory: %s', path, extra={'event': 'input.missing', 'error_code': 'INPUT_NOT_FOUND', 'action': 'Verify the selected ligand directory.'})
           raise FileNotFoundError(f'Ligand directory does not exist: {path}')

       self.ligandpath = path


    def set_receptorpath(self, path) -> None:
       if not os.path.exists(directory(path)):
           self.logger.error('Missing receptor directory: %s', path, extra={'event': 'input.missing', 'error_code': 'INPUT_NOT_FOUND', 'action': 'Verify the selected receptor directory.'})
           raise FileNotFoundError(f'Receptor directory does not exist: {path}')

       self.receptorpath = path


    def set_complexpath(self, path) -> None:
       if not os.path.exists(directory(path)):
           self.logger.error('Missing complex directory: %s', path, extra={'event': 'input.missing', 'error_code': 'INPUT_NOT_FOUND', 'action': 'Verify the selected complex directory.'})
           raise FileNotFoundError(f'Complex directory does not exist: {path}')

       self.complexpath = path


    def set_outputpath(self, path) -> None:
        if not os.path.exists(directory(path)):
            os.makedirs(directory(path), exist_ok=True)

        if not os.path.exists(directory(path)):
            raise OSError(f'Output directory could not be created: {path}')

        self.outputpath = path



    def generate_docking_script(self, input_template:str, output_script:str, **kwargs):
        """
        A generic function to generate docking scripts by replacing placeholders in the
        input template with provided keyword arguments. Such a function is useful to prepare
        different types of docking scripts for different docking softwares.

        Args:
            input_template (str): The path to the input template file.
            output_script (str): The path to the output script file.
            **kwargs: Keyword arguments used to replace placeholders in the input template.

        Raises:
            Exception: If any error occurs during the generation of the docking script, a log file will
            be created and posted in the log folder.

        """
        try:
            with open(resolve_path(input_template), 'r') as template_file:
                template_content = template_file.read()

            preparation = kwargs.pop('_preparation', None)
            record = kwargs.pop('_record', None)
            if preparation is not None:
                from biomolexplorer.redocking_config import configure_template
                template_content = configure_template(template_content, Path(input_template).name, preparation, record)
            from biomolexplorer.templates import quiet_template
            template_content=quiet_template('vina/config.template' if input_template.endswith('vina/config.template') else input_template,template_content)
            config_content = template_content.format(**kwargs)
            if 'runscript biomol_ligand_cleanup.py' in config_content:
                (Path(output_script).parent/'biomol_ligand_cleanup.py').write_text(
                    'import chimera\nimport DockPrep\n'
                    'DockPrep.prep(chimera.openModels.list(modelTypes=[chimera.Molecule]), '
                    'addHFunc=None, addCharges=False, delSolvent=False, delIons=False, '
                    'delLigands=False, delAltLocs=True, incompleteSideChains=False, '
                    'mutateMSE=False, mutate5BU=False, mutateUMS=False, mutateCSL=False, nogui=True)\n')
            if 'runscript biomol_receptor_cleanup.py' in config_content:
                helper=Path(output_script).parent/'biomol_receptor_cleanup.py'
                helper.write_text('import chimera\nimport DockPrep\n'
                    'from chimera.resCode import protein3to1\n'
                    'for model in chimera.openModels.list(modelTypes=[chimera.Molecule]):\n'
                    '    for residue in list(model.residues):\n'
                    '        if not residue.isHet or residue.type not in protein3to1 or residue.type in arguments: continue\n'
                    '        linked = any(n.residue != residue and n.residue.type in protein3to1 for a in residue.atoms for n in a.neighbors)\n'
                    '        if not linked:\n'
                    '            chimera.replyobj.info("Excluded non-polymer amino-acid additive: %s\\n" % residue.oslIdent())\n'
                    '            model.deleteResidue(residue)\n'
                    'DockPrep.prep(chimera.openModels.list(modelTypes=[chimera.Molecule]), '
                    'addHFunc=None, addCharges=False, delSolvent=False, delIons=False, '
                    'delLigands=False, delAltLocs=True, nogui=True)\n')
            # Classic Chimera's open/write parsers consume the remaining path
            # verbatim, including spaces. Shell quoting becomes part of its filename.

            with open(output_script, 'w') as output_file:
                output_file.write(config_content)

        except Exception as e:
            self.logger.error(f'during to perform {input_template} -> {output_script} in generate_docking_script function', exc_info=True)
            raise

        finally:

            sync_directory(directory(self.outputpath))




    def perform_subprocess(self, command, local_path=None, shell=False, check=True):
        """Execute a scientific tool; errors propagate to the job boundary.

        The shell argument remains accepted for legacy callers; execution always
        uses argv and explicit showbox input redirection.
        """
        return run_command(command, cwd=local_path, check=check)


    def prepare_on_chimera(self, filename:str) -> bool:
        """
        Some docking steps need to perform a specific input data file extension, and some features are
        required to report that. So, this function prepares such input data files when necessary,
        using the chimera shell execution to perform ones.

        Args:
            filename (str): The input file to be prepared for docking.

        Raises:
            Exception: If any error occurs during the preparation of the input file, a log file will
            be created and posted in the log folder.

        """
        try:
            import re
            script = Path(directory(self.outputpath), filename)
            entry=getattr(self,'_chimera_ligand_settings',{}).get(filename)
            if entry is not None and '\naddh\n' in script.read_text():
                from biomolexplorer.chimera_charge_repair import local_definition,charge_commands
                record,settings=entry
                if any(local_definition(self.complexpath,record[1],settings).values()):
                    source=script.read_text()
                    repairs=[line for line in charge_commands(self.complexpath,record,settings) if line not in source]
                    if repairs:script.write_text(source.replace('\naddh\n','\n'+'\n'.join(repairs)+'\naddh\n',1))
            outputs = re.findall(r'^write format (pdb|mol2) #0 (.+)$', script.read_text(), re.M)
            command = ['chimera', '--nogui', filename]
            with log_context(pair=Path(filename).stem.removeprefix('prepare_ligand_').removeprefix('prepare_receptor_').removeprefix('prepare_complex_')):
                try:
                    self.perform_subprocess(command, self.outputpath)
                except RuntimeError as error:
                    entry=getattr(self,'_chimera_ligand_settings',{}).get(filename)
                    if entry is None or 'Failure running ANTECHAMBER' not in str(error):raise
                    from biomolexplorer.chimera_charge_repair import charge_commands
                    record,settings=entry
                    repairs=charge_commands(self.complexpath,record,settings)
                    source=script.read_text()
                    if not repairs or '\naddh\n' not in source or all(line in source for line in repairs):raise
                    script.write_text(source.replace('\naddh\n','\n'+'\n'.join(repairs)+'\naddh\n',1))
                    event(self.logger,'ligand.typing.repaired','Retrying Chimera with CCD-confirmed atom types.',pair=filename)
                    self.perform_subprocess(command,self.outputpath)
            for file_format, output in outputs:
                self.validate_prepared_file(Path(directory(self.outputpath), output.strip()), file_format)
            script.unlink(missing_ok=True)
            return True

        except Exception as e:
            self.logger.error('Chimera preparation failed: %s', e, exc_info=True,
                extra={'event': 'preparation.failed', 'configuration': filename, **diagnose_exception(e)})
            raise

        finally:

            sync_directory(directory(self.outputpath))

    @staticmethod
    def validate_prepared_file(path, file_format):
        path = Path(path)
        if not path.is_file() or not path.stat().st_size:
            error = ValueError(f'Arquivo preparado ausente ou vazio: {path.name}. Consulte o log da ferramenta.')
            error.error_code = 'PREPARED_OUTPUT_INVALID'
            error.action = 'Inspect the preparation tool output and retained script; verify that the expected output contains atoms.'
            raise error
        text = path.read_text()
        if file_format in ('pdb', 'pdbqt'):
            atoms = any(line.startswith(('ATOM  ', 'HETATM')) for line in text.splitlines())
        elif file_format == 'mol2':
            atoms = '@<TRIPOS>ATOM' in text and bool(text.split('@<TRIPOS>ATOM', 1)[1].split('@<TRIPOS>', 1)[0].strip())
        else:
            return
        if not atoms:
            error = ValueError(f'Arquivo preparado sem átomos: {path.name}. Consulte o log da ferramenta.')
            error.error_code = 'PREPARED_OUTPUT_INVALID'
            error.action = 'Inspect the preparation tool output and retained script; verify that the expected output contains atoms.'
            raise error





    def prepare_on_obabel(self, inputfile:str, outputfile:str, params:Optional[List[Tuple[str, str]]]=[], input_format:Optional[str]='mol2', output_format:Optional[str]='pdb') -> bool:
        """
        Some docking steps need to perform a specific input data file extension, and some features are
        required to report that. So, this function prepares such input data files when necessary,
        using the obabel shell execution to perform ones.

        Args:
            inputfile (str): The input file to be prepared for docking.
            outputfile (str): The output file to be prepared for docking.
            params (List[Tuple[str, str]]): The list of parameters to be used in the obabel command.
            input_format (str): The input file format.
            output_format (str): The output file format.

        Raises:
            Exception: If any error occurs during the preparation of the input file, a log file will
            be created and posted in the log folder.

        """
        try:

            input_file = directory(self.outputpath) + inputfile
            if input_format=='mol2':
                from biomolexplorer.chimera_mol2 import normalize_carboxylates
                normalize_carboxylates(input_file)
                if inputfile.endswith('.lig.mol2'):
                    entry=next((value for name,value in getattr(self,'_chimera_ligand_settings',{}).items()
                        if name.removeprefix('prepare_ligand_').removesuffix('.com')==inputfile.removesuffix('.lig.mol2')),None)
                    if entry:
                        from biomolexplorer.dock6_validation import validate_mol2
                        from biomolexplorer.chimera_charge_repair import local_definition
                        if any(local_definition(self.complexpath,entry[0][1],entry[1]).values()):
                            from biomolexplorer.chimera_mol2 import canonicalize_component
                            canonicalize_component(input_file,self.complexpath,*entry)
                        else:
                            try:validate_mol2(input_file)
                            except ValueError:
                                from biomolexplorer.chimera_mol2 import canonicalize_component
                                canonicalize_component(input_file,self.complexpath,*entry)
                        validate_mol2(input_file)
            if input_format != 'smi':
                self.validate_prepared_file(input_file, input_format)

            if input_format != 'smi':
                command = ['obabel', '-i', input_format, input_file, '-o', output_format, '-O', outputfile]
            else:
                command = ['obabel', f'-:{inputfile}', '-o', output_format, '-O', outputfile, '--gen3D']


            for param, value in params:
                command.append(f'-{param}')
                if value:
                    command.append(str(value))

            self.perform_subprocess(command, self.outputpath)
            self.validate_prepared_file(Path(directory(self.outputpath), outputfile), output_format)
            return True

        except Exception as e:
            self.logger.error(f'during to perform {inputfile} to {outputfile} converter in prepare_on_obabel function', exc_info=True)
            raise

        finally:

            sync_directory(directory(self.outputpath))





    def calculate_ligand_centerofmass(self, inputfile:str, ligand:str):
        """
        This function calculates the center of mass of a given ligand using the PyMOL software.

        Args:
            inputfile (str): The input file containing the ligand.
            ligand (str): The name of the ligand in the input file.

        Returns:
            List: The center of mass of the ligand.

        Raises:
            Exception: If any error occurs during the calculation of the ligand center of mass, a log file
            will be created and posted in the log folder.
        """
        try:

            file = directory(self.outputpath) + inputfile
            cmd.reinitialize()
            cmd.load(file, 'ligand')
            cmd.select('ligante', 'resn ' + ligand)
            return cmd.centerofmass('ligante')

        except Exception as e:
            self.logger.error(f'during to perform {inputfile} -> {ligand} ligand in calculate_ligand_centerofmass function', exc_info=True)
            return None



    def retrieve_centerofmass_dataset(self, pathname:str, receptor:str, ligand:str, resnum:str, chain:str) -> List:
        """
        This function retrieves the center of mass of a given ligand from a dataset.

        Args:
            pathname (str): The path to the dataset file.
            receptor (str): The receptor name.
            ligand (str): The ligand name.
            resnum(str): The residue number.
            chain(str): The chain identifier.

        Returns:
            List: The center of mass of the ligand.

        Raises:
            Exception: If any error occurs during the retrieval of the ligand center of mass, a log file
        """
        try:
            if self.centers is None:
                f1   = fileHandling(input_path=pathname, output_path=self.outputpath)
                self.centers = f1.csv_to_dataframe('centers')

            center_idx = f'{receptor}_{ligand}_{resnum}{chain}'
            legacy_idx = f'{receptor}_{ligand}_{resnum}_{chain}'
            column=center_idx if center_idx in self.centers else legacy_idx
            return self.centers[column].values.tolist() if column in self.centers else None

        except Exception as e:
            self.logger.error(f'during to perform {receptor} -> {ligand }in retrieve_centerofmass_dataset function', exc_info=True)
            raise



    def process_in_parallel(self, method_name:str, args_list:list, process_by_threads:Optional[bool]=False):
        """
        A broad function to perform a parallel processing of a given method with a list of arguments. In this function, it is possible
        to choose between processing by threads (True) or by processes (False), modifying the process_by_threads input variable before execution.

        Args:
            method_name (str): The method name to be processed in parallel.
            args_list (list): The list of arguments to be processed in parallel.
            nprocess (int): The number of processes to be used in parallel.
            process_by_threads (bool): If True, the processing will be done using threads.

        Returns:
            List: The results of the parallel processing.

        Raises:
            Exception: If any error occurs during the parallel processing, a log file will be created
            and posted in the log folder.
        """


        method   = getattr(self, method_name)

        if process_by_threads:
            with ThreadPoolExecutor(max_workers=self.nprocess) as executor:
                futures = [executor.submit(method, *args) for args in args_list]
                results = [future.result() for future in futures]

        else:
            with ProcessPoolExecutor(max_workers=self.nprocess) as executor:
                futures = [executor.submit(method, *args) for args in args_list]
                results = [future.result() for future in futures]

        return results



    def prepare_for_docking(self, pdb_codes:list, charge_type:str, pH:float, redefine_centerofmass:bool, preparation_pairs=None, _prepared_receptors=None) -> bool:
        """
        This function prepares the input files for docking using Chimera and Open Babel software. The input files are prepared in a specific
        format required to perform docking with different docking software. Templates for the input files are available in the src/biomolexplorer/resources folder.
        Some limitations are present in the preparation of the input files, such as the need to preserve declared variables because they are used
        to replace placeholders in the output scripts prepared.

        Args:
            pdb_codes (List[Tuple[str, str]]): The list of PDB codes to be used in the docking preparation. If None is provided, such codes are
            extracted utilizing the default PDB pathway provided in the complexpath variable to retrieve the PDB codes.
            pH (float): The pH value to be used for protonation. If It is not explicitly described, the default value used is 7.4.
            redefine_centerofmass (bool): If True, the function will calculate the center of mass of the ligand.

        Raises:
            Exception: If any error occurs during the preparation of the input files for docking, a log file will be created and posted in the log folder.

        """
        try:

            from biomolexplorer.redocking_config import pair_key, DEFAULTS
            from biomolexplorer.redocking_config import chimera_pairs
            preparation_pairs=chimera_pairs(pdb_codes,preparation_pairs)
            event(self.logger, 'preparation.started', 'Preparing selected receptor/ligand complexes')
            centers   = {}
            complexes = []
            receptors = []
            ligands   = []
            self._chimera_ligand_settings={}


            for pdb in pdb_codes:
                receptor = pdb[0]
                ligand   = pdb[1]
                resnum   = pdb[2]
                chain    = pdb[3]
                preparation = (preparation_pairs or {}).get(pair_key(pdb), {})
                preparation = dict(preparation, ligand={'charge_type': charge_type, **preparation.get('ligand', {})})
                chain_id = chain
                complex_preparation = None
                if preparation_pairs is not None:
                    complex_preparation = preparation

                if (_prepared_receptors is None or (receptor,chain) not in _prepared_receptors) and f'prepare_complex_{receptor}_{chain}.com' not in complexes:
                    self.generate_docking_script(input_template='src/scripts/chimera/prepare_complex.template',
                                                output_script=directory(self.outputpath) + f'prepare_complex_{receptor}_{chain}.com',
                                                pdb_code=directory(self.complexpath) + receptor,
                                                input_complex=directory(self.complexpath) + receptor,
                                                chain=chain_id,
                                                output_complex=directory(self.outputpath) + f'{receptor}_{chain}', _preparation=complex_preparation, _record=pdb)
                    complexes.append(f'prepare_complex_{receptor}_{chain}.com')


                if (_prepared_receptors is None or (receptor,chain) not in _prepared_receptors) and f'prepare_receptor_{receptor}_{chain}.com' not in receptors:
                    self.generate_docking_script(input_template='src/scripts/chimera/prepare_receptor.template',
                                                output_script=directory(self.outputpath) + f'prepare_receptor_{receptor}_{chain}.com',
                                                input_complex=directory(self.outputpath) + f'{receptor}_{chain}' + '.complex',
                                                receptor=directory(self.outputpath) + f'{receptor}_{chain}', _preparation=preparation if preparation_pairs is not None else None)
                    receptors.append(f'prepare_receptor_{receptor}_{chain}.com')


                extention = f'{receptor}_{ligand}_{resnum}{chain}'
                self.generate_docking_script(input_template='src/scripts/chimera/prepare_ligand.template',
                                             output_script=directory(self.outputpath) + f'prepare_ligand_{extention}.com',
                                             input_complex=directory(self.outputpath) + f'{receptor}_{chain}' + '.complex',
                                             resnum=resnum,
                                             chain=chain,
                                             charge_type=charge_type,
                                             input_ligand=directory(self.outputpath) + f'{extention}',
                                             output_ligand=directory(self.outputpath) + f'{extention}', _preparation=preparation if preparation_pairs is not None else None)
                ligands.append(f'prepare_ligand_{extention}.com')
                self._chimera_ligand_settings[f'prepare_ligand_{extention}.com']=(pdb,preparation)



            args = [(file,) for file in complexes]
            if _prepared_receptors is None:self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)
            else:
                for (file,) in args:self.prepare_on_chimera(file)

            args = [(file,) for file in receptors]
            if _prepared_receptors is None:self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)
            else:
                for (file,) in args:self.prepare_on_chimera(file)

            converted_receptors = set()
            for pdb in pdb_codes:
                chain = pdb[3]
                receptor_key = (pdb[0], chain)
                if _prepared_receptors is not None and receptor_key in _prepared_receptors:continue
                if receptor_key in converted_receptors:
                    continue
                options = dict(DEFAULTS, **(preparation_pairs or {}).get(pair_key(pdb), {}).get('receptor', {}))
                self.prepare_on_obabel(f'{pdb[0]}_{chain}.dockprep.mol2', f'{pdb[0]}_{chain}.dockprep.pdbqt', ([('p',str(pH))] if options['add_hydrogens'] else []) + [('xr','')], input_format='mol2', output_format='pdbqt')
                converted_receptors.add(receptor_key)
                if _prepared_receptors is not None:_prepared_receptors.add(receptor_key)
                checkpoint=getattr(self,'receptor_checkpoint',None)
                if checkpoint:checkpoint(pdb)

            args = [(file,) for file in ligands]
            if _prepared_receptors is None:self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)
            else:
                for (file,) in args:self.prepare_on_chimera(file)

            tmp_codes = []
            for pdb in pdb_codes:
                chain = pdb[3]
                key = f'{pdb[0]}_{pdb[1]}_{pdb[2]}{chain}'
                options = dict(DEFAULTS, **(preparation_pairs or {}).get(pair_key(pdb), {}).get('ligand', {}))
                self.prepare_on_obabel(f'{key}.lig.mol2', f'{key}.lig.pdbqt', [('p',str(pH))] if options['add_hydrogens'] else [], input_format='mol2', output_format='pdbqt')
                if redefine_centerofmass:
                    center = self.calculate_ligand_centerofmass(f'{key}.lig.pdb', pdb[1])
                    if center is None or len(center)!=3 or not all(math.isfinite(float(v)) for v in center):
                        raise ValueError(f'Não foi possível calcular o centro do ligante {key}. Verifique o arquivo preparado e o log da etapa.')
                    centers[key] = center
                    tmp_codes.append(tuple(pdb))


            if redefine_centerofmass:
                f1 = fileHandling(input_path=self.outputpath, output_path=self.outputpath)
                tmp = f1.csv_to_dataframe('centers').to_dict(orient='list') if f1.isFile('centers')[0] else {}
                tmp.update(centers)
                df = DataFrame(tmp)
                f1.dataframe_to_csv('centers', df)

            event(self.logger, 'preparation.succeeded', 'Complex preparation completed', artifacts=len(tmp_codes))
            return tmp_codes if redefine_centerofmass else list(pdb_codes)


        except Exception:
            self.logger.exception('Falha na preparação dos complexos para docking.')
            raise

        finally:
            sync_directory(directory(self.outputpath))

    def recover_better_conforms_of_vina(self, charge_type:str, filename:Optional[list]=None, molecules_dataset:Optional[str]=None, preparation_pairs=None):
        """
        This function recover better conformations for each ligand performed by AutoDock Vina software. The input files are prepared in a specific
        format required to perform docking with different docking software. Templates for the input files are available in the src/biomolexplorer/resources folder.
        Some limitations are present in the preparation of the input files, such as the need to preserve declared variables because they are used
        to replace placeholders in the output scripts prepared.

        Args:
            filename (list): The list of files to be prepared for better conformation recovery. If None is provided, such files are extracted
            utilizing the default ligand pathway provided in the ligandpath variable to retrieve the ligand files.
            molecules_dataset (str): The path to the dataset file containing the molecules to be used in the docking analysis.

        Raises:
            Exception: If any error occurs during the recovery of better conformations, a log file will be created and posted in the log folder.

        """
        try:

            files = [f for f in os.listdir(directory(self.ligandpath)) if f.endswith('.pdbqt')] if filename == None else filename
            conv      = MolConverter(input_path=self.ligandpath, output_path=self.outputpath)
            explorer  = MolExplorer(input_path=self.ligandpath, output_path=self.outputpath)

            for input_file in files:

                pdb_str_content = conv.extract_pdb_to_pdbqt(input_file, start_index='MODEL 1', end_index='MODEL 2', pdb_filename=input_file.rsplit('.')[0]+'.pdb')


                if not explorer.is_fragmented(pdb_str_content):

                    self.generate_docking_script(input_template='src/scripts/chimera/prepare_better_conform.template',
                                             output_script=directory(self.outputpath) + f'prepare_better_conform_{input_file}.com',
                                             ligand=input_file.rsplit(".")[0], charge_type=charge_type,
                                             _preparation=next((config for key, config in (preparation_pairs or {}).items() if input_file.startswith('_'.join(key.split('|')[:3]) + key.split('|')[3] + '.')), None))


            args = [(f'prepare_better_conform_{input_file}.com',) for input_file in files]
            self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)


            files = [f for f in os.listdir(directory(self.outputpath)) if f.endswith('.lig.mol2')]
            for file in files:
                with open(directory(self.outputpath) + file, 'r') as fp:
                    lines = ''.join(fp.readlines())
                    if lines.find('nan') >= 0:
                        os.remove(directory(self.outputpath) + file)


            if molecules_dataset:
                fx = fileHandling(input_path=molecules_dataset, output_path=molecules_dataset)
                files = [f.replace('.lig.mol2','').split('_')[1] for f in os.listdir(directory(self.outputpath)) if f.endswith('.lig.mol2')]
                df = fx.csv_to_dataframe(self.mol_filename)
                data = df[df['molecule_chembl_id'].isin(files)]
                fx.dataframe_to_csv(self.mol_filename, data)
                data = set(files) - set(data['molecule_chembl_id'].tolist())
                path = directory(self.outputpath)
                files = [os.remove(path + f) for f in os.listdir(path) if f.replace('.lig.mol2','').split('_')[1] in data]


        except Exception as e:
            self.logger.error(f'during to perform {filename} in recover_better_conformation function', exc_info=True)
            raise

        finally:
            sync_directory(directory(self.outputpath))

            pdbs = [f for f in os.listdir(directory(self.outputpath)) if f.endswith('.pdb')]
            [os.remove(directory(self.outputpath) + file) for file in pdbs]

            mols = [f for f in os.listdir(directory(self.outputpath)) if f.endswith('(2).lig.mol2')]
            [os.remove(directory(self.outputpath) + file) for file in mols]

            sync_directory(directory(self.outputpath))

class DockVina(Docking):


    def __init__(self, ligand_input_path:Optional[str]=None, receptor_input_path:Optional[str]=None,
                 complex_input_path:Optional[str]=None, output_path:Optional[str]=None, mol_filename:Optional[str]='molecules',
                 pdb_codes:Optional[Tuple[str, str, str, str]]=None, centerofmasspath:Optional[str]=None,
                 sizeof_box:Optional[List]=[24,24,24], exhaustiveness:Optional[int]=20, num_modes:Optional[int]=10) -> None:


        super().__init__(ligand_input_path, receptor_input_path, complex_input_path, output_path, mol_filename=mol_filename)

        self.__pdb_codes        = pdb_codes
        self.__centerofmasspath = centerofmasspath
        self.__sizeof_box       = sizeof_box
        self.__exhaustiveness   = exhaustiveness
        self.__num_modes        = num_modes


    def prepare_compounds_for_vina(self, pH:Optional[float]=7.4):
        """
        This function prepares the input files for docking using the AutoDock Vina. The input files are prepared in a specific format (i.e. PDBQT) required
        to perform docking. For that, the molecules.csv file, available in the graph/data path, retrieves the molecules in SMILES format to be prepared for
        docking. If necessary, the molecules file can be updated to introduce new molecules to be prepared for docking. For that, include the new lines in
        the molecules.csv file according to the following format: molecule_chembl_id, canonical_smiles. In moleculle_chmbl_id, the molecule identifier is
        defined as an investigated molecule, so you may want to use a specific/particular code to represent the newly introduced lines in the file.

        Args:
            pH (float): The pH value to be used for protonation. If It is not explicitly described, the default value used is 7.4.

        Raises:
            Exception: If any error occurs during the preparation of the input files for docking, a log file will be created and posted in the log folder.

        """
        df = None
        f1 = None
        try:

            f1 = fileHandling(input_path=self.ligandpath, output_path=self.outputpath)
            df = f1.csv_to_dataframe(self.mol_filename)

            df['molecule_chembl_id'] = df['molecule_chembl_id'].astype(str)
            data = df[['molecule_chembl_id', 'canonical_smiles']].to_records(index=False)

            for chemblid, smiles in data:
                self.prepare_on_obabel(inputfile=smiles, outputfile=chemblid + '.lig.pdbqt', input_format="smi", output_format='pdbqt', params=[('p',str(pH))])


        except Exception as e:
            self.logger.error('Compound preparation failed', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            raise

        finally:
            sync_directory(directory(self.outputpath))

            mols = [f.split('.lig.pdbqt')[0] for f in os.listdir(directory(self.outputpath)) if f.endswith('.lig.pdbqt')]
            if df is not None and f1 is not None:
                df = df[df['molecule_chembl_id'].isin(mols)]
                f1.dataframe_to_csv(self.mol_filename, df)

            sync_directory(directory(self.outputpath))

    def perform_vina_evaluation(self):
        """
        This function performs the AutoDock Vina evaluation for each prepared ligand in the ligand path folder. The assessment
        is performed sequentially. The results are stored in the output path folder in PDBQT files named according to the
        identifications reported in molecules.csv.

        Raises:
            Exception: If any error occurs during the evaluation of the prepared ligands for docking, a log file will
            be created and posted in the log folder.

        """

        files_to_perform = [f for f in os.listdir(directory(self.outputpath)) if f.endswith('.vina')]

        from biomolexplorer.vina_batch import run_resilient
        return run_resilient([Path(self.outputpath)/name for name in files_to_perform],self.outputpath)



    def redocking(self, pH:float):
        """
        The redocking analysis is conducted using the AutoDock Vina software with default parameters. The search box size is set to [24,24,24], the exhaustiveness
        is set to 20, and the number of result modes is set to 10. The parameters can be modified based on user conditions, and results are stored in a series of PDBQT files,
        named according to the identifiers reported in molecules.csv and stored in the output pathway. During this stage, the Root Mean Square Deviation (RMSD) is calculated
        to evaluate the quality of the docking analysis. For RMSD calculation, a CSV file is generated in the output path folder, containing the RMSD values for each docking analysis.

        Raises:
            Exception: If any error occurs during the redocking analysis, a log file will be created and posted in the log folder.

        """
        try:

            desc = Descriptors()
            results = []
            pdb_codes = DataFrame(self.__pdb_codes, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN'])
            pdb_codes['RESNUM'] = pdb_codes['RESNUM'].astype(str)
            pdb_codes = pdb_codes.to_records(index=False)

            for idx, (receptor, ligand, resnum, chain) in enumerate(pdb_codes):
                composite = f'{receptor}_{ligand}_{resnum}{chain}'

                center = self.retrieve_centerofmass_dataset(self.ligandpath, receptor, ligand, resnum, chain)
                if center == None:
                    raise ValueError(f'Centro do ligante ausente ou inválido: {composite}')

                self.generate_docking_script(input_template='src/scripts/vina/config.template',
                                            output_script=directory(self.outputpath) + f'{composite}.vina',
                                            receptor=directory(self.receptorpath) + f'{receptor}_{chain}.dockprep.pdbqt',
                                            ligand=directory(self.ligandpath) + f'{composite}' + '.lig.pdbqt',
                                            center_x=center[0],
                                            center_y=center[1],
                                            center_z=center[2],
                                            size_x=self.__sizeof_box[0],
                                            size_y=self.__sizeof_box[1],
                                            size_z=self.__sizeof_box[2],
                                            out=f'{composite}.lig.pdbqt',
                                            exhaustiveness=self.__exhaustiveness,
                                            num_modes=self.__num_modes)



            self.perform_vina_evaluation()

            for receptor, ligand, resnum, chain in pdb_codes:
                composite = f'{receptor}_{ligand}_{resnum}{chain}'
                iligand  = directory(self.ligandpath) + f'{composite}' + '.lig.pdbqt'
                vina_model = directory(self.outputpath) + f'{composite}' + '.lig.pdbqt'
                if not os.path.isfile(iligand) or not os.path.isfile(vina_model):
                    raise ValueError(f'Saída de redocking ausente: {composite}. Consulte o log do Vina.')
                rmsd_value = desc.calcRMSD(iligand, vina_model)
                if not math.isfinite(float(rmsd_value)):
                    raise ValueError(f'RMSD inválido no redocking: {composite}')
                results.append((f'{receptor}', f'{ligand}', f'{resnum}', f'{chain}', rmsd_value))
                callback=getattr(self,'pair_checkpoint',None)
                if callback:callback([receptor,ligand,resnum,chain],rmsd_value)


            rmsd = DataFrame(results, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RMSD'])
            pdb_codes = DataFrame(self.__pdb_codes, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION'])
            pdb_codes['RESNUM'] = pdb_codes['RESNUM'].astype(str)
            pdb_codes = pdb_codes.merge(rmsd, on=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN'], how='left')

            f1   = fileHandling(output_path=self.complexpath)
            f1.dataframe_to_csv('pdb_codes', pdb_codes)


        except Exception as e:
            self.logger.error('during to perform the docking function', exc_info=True)
            raise



    def docking(self, base_selected_mols:str):

        """
        This function performs the docking analysis using the AutoDock Vina software. The analysis is conducted sequentially, considering the number of processes
        available in the molecules.csv. The search box size is set to [24,24,24], the exhaustiveness is set to 20, and the number of result modes is set to 10. The parameters
        can be modified based on user conditions, and results are stored in a series of PDBQT files, named according to the identifiers reported in molecules.csv.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.

        """

        try:

            molecules = [f.rsplit('.lig.pdbqt')[0] for f in os.listdir(directory(self.ligandpath)) if f.endswith('.lig.pdbqt')]
            selected=getattr(self,'selected_ligand_files',None)
            if selected is not None:
                allowed={Path(p).name.removesuffix('.lig.pdbqt') for p in selected}
                molecules=[code for code in molecules if code in allowed]
            configurations=[]

            self.__pdb_codes   = [(pdb[0], pdb[1], pdb[2], pdb[3]) for pdb in self.__pdb_codes]

            for receptor, ligand, resnum, chain in self.__pdb_codes:
                center    = self.retrieve_centerofmass_dataset(self.__centerofmasspath, receptor, ligand, resnum, chain)
                if center is None:raise ValueError(f'Docking center is missing for {receptor}_{ligand}_{resnum}{chain}.')
                reference=f'{receptor}_{ligand}_{resnum}{chain}'
                prefix=reference+'_'
                # File existence alone is not proof of a compatible experiment.
                # The per-calculation content index decides what can be reused.
                for mol in sorted(set(molecules)):
                     self.generate_docking_script(input_template='src/scripts/vina/config.template',
                                            output_script=directory(self.outputpath) + f'{reference}_{mol}.vina',
                                            receptor=directory(self.receptorpath) + f'{receptor}_{chain}' + '.dockprep.pdbqt',
                                            ligand=directory(self.ligandpath) + mol + '.lig.pdbqt',
                                            center_x=center[0],
                                            center_y=center[1],
                                            center_z=center[2],
                                            size_x=self.__sizeof_box[0],
                                            size_y=self.__sizeof_box[1],
                                            size_z=self.__sizeof_box[2],
                                            out=f'{reference}_{mol}.lig.pdbqt',
                                            exhaustiveness=self.__exhaustiveness,
                                            num_modes=self.__num_modes)

                     sync_directory(directory(self.outputpath))

                     configuration = Path(self.outputpath) / f'{reference}_{mol}.vina'
                     configurations.append(configuration)

            # The Vina Python extension may terminate the interpreter for a
            # malformed ligand.  A checkpointed child retains map reuse in the
            # common case and lets this worker resume after such native aborts.
            from biomolexplorer.vina_batch import run_resilient
            return run_resilient(configurations,self.outputpath)


        except Exception as e:
            self.logger.error('during to perform the docking function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            raise

        finally:

            sync_directory(directory(self.outputpath))

            [os.remove(directory(self.outputpath) + file) for file in os.listdir(directory(self.outputpath)) if file.endswith('.vina')]
            [os.remove(directory(self.outputpath) + file) for file in os.listdir(directory(self.outputpath)) if file.endswith('(2).pdbqt')]
            self.centers = None




class Dock6(Docking):


    def __init__(self, dock6_path:Optional[str]='', ligand_input_path:Optional[str]=None, receptor_input_path:Optional[str]=None,
                 base_output_path:Optional[str]=None, pdb_code:Optional[str]=None, density:Optional[float]=0.5,
                 radius:Optional[float]=1.4, distance:Optional[float]=10.0, max_residues:Optional[int]=50,
                 conformer_search_type:Optional[Literal['flex', 'rigid']] = 'flex', mol_filename:Optional[str]='molecules',
                 binding_site_center=None) -> None:

        # Legacy DOCK6 helpers truncate long paths and split whitespace.
        import tempfile
        import shutil
        self.runtime_path=tempfile.mkdtemp(prefix='bme-',dir='/tmp')
        self.export_path=base_output_path
        runtime=Path(self.runtime_path)
        shutil.copytree(ligand_input_path,runtime/'lig')
        shutil.copytree(receptor_input_path,runtime/'rec')
        from biomolexplorer.chimera_mol2 import normalize_carboxylates
        for folder in ('lig','rec'):
            for mol2 in (runtime/folder).glob('*.mol2'):normalize_carboxylates(mol2)
        super().__init__(ligand_input_path=str(runtime/'lig'),
                         receptor_input_path=str(runtime/'rec'),
                         output_path=str(runtime/'out'),
                         mol_filename=mol_filename)

        self.__dock6_path            = directory(dock6_path)
        self.__base_output_path      = str(runtime/'out')
        self.__pdb_code              = pdb_code
        self.__density               = density
        self.__radius                = radius
        self.__distance              = distance
        self.__max_residues          = max_residues
        self.__conformer_search_type = conformer_search_type
        self.binding_site_center = binding_site_center

    def export_results(self):
        import shutil
        shutil.copytree(self.__base_output_path,self.export_path,dirs_exist_ok=True)
        shutil.rmtree(self.runtime_path)

    def set_compound_report(self, dataset, report):
        from biomolexplorer.docking_data import read_compounds
        self.compound_dataset=Path(dataset)
        self.compound_rows={row['molecule_chembl_id']:(index,row) for index,row in enumerate(read_compounds(dataset),2)}
        self.failure_report=report
        self.excluded_compounds=set()

    def _ensure_compounds(self):
        if not hasattr(self,'failure_report'):
            from biomolexplorer.preparation_failures import CURRENT_REPORT,PreparationFailures
            self.failure_report=CURRENT_REPORT.get() or PreparationFailures()
            self.compound_dataset=Path(self.ligandpath)/'compounds.csv'
            self.compound_rows={}
            self.excluded_compounds=set()

    def _compound_guard(self, code, step):
        self._ensure_compounds()
        index,row=self.compound_rows.get(code,(None,{'molecule_chembl_id':code}))
        names={code+'.lig.min_scored.mol2',code+'_scored.mol2',code+'_footprint_scored.txt',code+'.pdf'}
        cleanup=[p for p in Path(self.__base_output_path).rglob('*') if p.is_file() and p.name in names]
        return self.failure_report.compound(self.compound_dataset,index,row,self.export_path,cleanup,step=step)

    def _run_compound_jobs(self, method, jobs, step):
        """Collect failures in the parent thread; never lose child-process reports."""
        from concurrent.futures import as_completed
        from biomolexplorer.dock6_validation import validate_mol2
        from biomolexplorer.footprints import footprint_rows
        from biomolexplorer.progress import report_progress
        self._ensure_compounds()
        checkpoints=self._checkpoints() if hasattr(self,'runtime_path') else None
        pending=[];keys={}
        for code,args,output in jobs:
            if checkpoints:
                key=checkpoints.compound_key(method,code);keys[code]=key
                def validate(role,path):
                    if output.suffix=='.mol2':validate_mol2(path,scored=True)
                    else:footprint_rows(path)
                if checkpoints.cache.restore(key,{'result':output},validate) is not None:continue
            pending.append((code,args,output))
        def execute_compound(code,args,output):
            if not checkpoints:return getattr(self,method)(*args)
            def validate(role,path):
                if output.suffix=='.mol2':validate_mol2(path,scored=True)
                else:footprint_rows(path)
            with checkpoints.cache.locked(keys[code]):
                if checkpoints.cache.restore(keys[code],{'result':output},validate) is not None:return
                getattr(self,method)(*args)
                validate('result',output)
                checkpoints.cache.remember(keys[code],{'result':output},dict(compound=code,step=step,
                    receptor_identifier=checkpoints.receptor_identifier))
        with ThreadPoolExecutor(max_workers=self.nprocess) as executor:
            futures={executor.submit(execute_compound,code,args,output):(code,output) for code,args,output in pending}
            for completed,future in enumerate(as_completed(futures),1):
                code,output=futures[future];success=False
                with self._compound_guard(code,step):
                    future.result()
                    if output.suffix=='.mol2':validate_mol2(output,scored=True)
                    else:footprint_rows(output)
                    success=True
                if not success:self.excluded_compounds.add(code)
                report_progress(step+': '+code,completed,len(jobs))
        self.failure_report.write(self.export_path)

    def _checkpoints(self):
        if not hasattr(self,'calculation_checkpoints'):
            from biomolexplorer.dock6_cache import Dock6Checkpoints
            self.calculation_checkpoints=Dock6Checkpoints(self,self.__dock6_path,self.__pdb_code,
                dict(density=self.__density,radius=self.__radius,distance=self.__distance,
                     center=self.binding_site_center),self.__base_output_path)
        return self.calculation_checkpoints

    def prepare_site(self):
        checkpoints=self._checkpoints()
        with checkpoints.cache.locked(checkpoints.site_key()):
            if checkpoints.restore_site():return
            self.prepare_surface();self.prepare_showbox();self.prepare_gridbox()
            checkpoints.remember_site()

    def validate_setup(self):
        from biomolexplorer.processes import ScientificToolError
        root=Path(self.__dock6_path)
        for name in ('dock6','sphgen','sphere_selector','showbox','grid'):
            file=root/'bin'/name
            if not file.is_file() or not os.access(file,os.X_OK):
                raise ScientificToolError('Executável DOCK6 ausente ou não executável: '+str(file),'TOOL_NOT_FOUND','Revise o diretório de instalação do DOCK6.')
        for name in ('vdw_AMBER_parm99.defn','flex.defn','flex_drive.tbl'):
            file=root/'parameters'/name
            if not file.is_file() or not file.stat().st_size:raise ValueError('Parâmetro global DOCK6 ausente ou vazio: '+str(file))
        self._ensure_compounds()
        from biomolexplorer.dock6_validation import compound_code,validate_mol2
        for file in sorted(Path(self.ligandpath).glob('*.mol2')):
            code=compound_code(file);success=False
            with self._compound_guard(code,'Validação do ligante MOL2 para DOCK6'):
                validate_mol2(file)
                success=True
            if not success:self.excluded_compounds.add(code)
        self.failure_report.write(self.export_path)

    def perform_subprocess(self, command, local_path=None, shell=False, check=True):
        import shlex
        arguments=shlex.split(command) if isinstance(command,str) else list(command)
        if arguments[0] in ('dock6','sphgen','sphere_selector','showbox','grid'):
            arguments[0]=str(Path(self.__dock6_path)/'bin'/arguments[0])
        return super().perform_subprocess(arguments,local_path,shell,check)

    def _grid_prefix(self, code):
        name='binding_site' if getattr(self,'binding_site_center',None) is not None else code
        return str(Path(self.__base_output_path)/'gridbox'/name)

    def _ligand_files(self, folder):
        from biomolexplorer.dock6_validation import compound_code
        self._ensure_compounds()
        return [p.name for p in sorted(Path(folder).glob('*.mol2')) if compound_code(p) not in self.excluded_compounds]



    def prepare_surface(self) -> None:

        """
        The following function outlines the steps to prepare the surface of a pdb_code using the DOCK 6 software.
        As default, the density is set to 0.5 and the radius to 1.4, but, if necessary, these parameters can be modified.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.

        """

        try:

            self.set_outputpath(f'{self.__base_output_path}/surface/' )

            os.remove(directory(self.outputpath) + self.__pdb_code +".sph") if os.path.exists(directory(self.outputpath) + self.__pdb_code +".sph") else None

            input = directory(self.receptorpath) + self.__pdb_code
            output = directory(self.outputpath) + self.__pdb_code

            summary = generate_surface(f'{input}.noH.pdb', f'{output}.dms',
                                       density=self.__density, probe_radius=self.__radius)
            self.logger.info('Native molecular surface: %s', summary)

            self.generate_docking_script(input_template='src/scripts/dock6/INSPH.template',
                                         output_script=directory(self.outputpath)+'INSPH',
                                         receptor=self.__pdb_code)

            command = f'sphgen -i INSPH -o OUTSPH'
            self.perform_subprocess(command, self.outputpath)

            #prepare the surface selectors
            self.set_outputpath(f'{self.__base_output_path}/surface/Molecules/')
            selector = directory(f'{self.__base_output_path}/surface/') + self.__pdb_code + '.sph'

            site=None
            if getattr(self,'binding_site_center',None) is not None:
                x,y,z=self.binding_site_center
                site=Path(self.outputpath)/'binding_site.mol2'
                site.write_text('@<TRIPOS>MOLECULE\nsite\n1 0 1 0 0\nSMALL\nNO_CHARGES\n'
                    f'@<TRIPOS>ATOM\n1 C {x} {y} {z} C.3 1 SITE 0.0\n'
                    '@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n1 SITE 1\n')

            files = self._ligand_files(self.ligandpath)
            if site:files=files[:1]
            for ligand in files:
                command = ['sphere_selector', selector, str(site) if site else directory(self.ligandpath) + ligand, str(self.__distance)]
                self.perform_subprocess(command, self.outputpath)
                from biomolexplorer.dock6_validation import compound_code
                selected=Path(self.outputpath)/'selected_spheres.sph'
                if not selected.is_file() or not selected.stat().st_size:raise ValueError('A seleção do sítio DOCK6 não produziu esferas.')
                if site:
                    spheres=[]
                    for line in selected.read_text().splitlines():
                        fields=line.split()
                        try:
                            if len(fields)>=5 and int(fields[0])>0 and all(math.isfinite(float(v)) for v in fields[1:5]) and float(fields[4])>0:spheres.append(fields)
                        except ValueError:pass
                    if not spheres:raise ValueError('Nenhuma esfera válida foi selecionada para o sítio DOCK6. Revise centro e distância.')
                os.rename(selected,Path(self.outputpath)/('binding_site.sph' if site else compound_code(ligand)+'.sph'))

        except Exception as e:
            self.logger.error(f'during to perform into prepare_surface function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            raise

        finally:

            output = f'{directory(self.__base_output_path)}/surface/'

            os.remove(output + "temp1.ms") if os.path.exists(output + "temp1.ms") else None
            os.remove(output + "temp2.sph") if os.path.exists(output + "temp2.sph") else None
            os.remove(output + "temp3.atc") if os.path.exists(output + "temp3.atc") else None

            sync_directory(output)

    def prepare_showbox(self):
        """
        This function details the steps to prepare the showbox of a receptor using the DOCK 6 software. The analysis is conducted
        sequentially, taking into account the number of processes specified in molecules.csv. The search box size is set to [24, 24, 24],
        the exhaustiveness is set to 20, and the number of result modes is set to 10. These parameters can be modified based on user
        requirements. The results are stored in a series of PDBQT files, named according to the identifiers in molecules.csv.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """
        try:

            self.set_outputpath(f'{self.__base_output_path}/showbox/')
            surface_inputpath = f'{self.__base_output_path}/surface/Molecules/'

            files = [f for f in os.listdir(directory(surface_inputpath)) if f.endswith('.sph')]
            for ligand in files:

                self.generate_docking_script(input_template='src/scripts/dock6/showbox.template',
                                             output_script=directory(self.outputpath) + f'{Path(ligand).stem}.in',
                                             in_surface='../surface/Molecules/' + ligand,
                                             out_surface=Path(ligand).stem+'.box.pdb')


            for ligand in files:
                command = ['showbox','<',Path(ligand).stem+'.in']
                self.perform_subprocess(command, self.outputpath)
                box=Path(self.outputpath)/(Path(ligand).stem+'.box.pdb')
                if not box.is_file() or not box.stat().st_size:raise ValueError('A etapa showbox não produziu a caixa do sítio: '+box.name)


        except Exception as e:
            self.logger.error(f'during to perform the prepare_showbox function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            raise

        finally:


            sync_directory(directory(self.outputpath))





    def perform_parallel_gridbox(self, showbox:str, tid:int):
        """
        This function outlines the steps to prepare the grid box of a receptor using the DOCK 6 software. The analysis is conducted
        in a parallel execution, taking into account the number of CPU's specified in the computational architecture.

        Args:
            receptor (str): The receptor file to perform the grid box.
            showbox (str): The showbox file to perform the grid box.
            tid (int): The identifier of the thread.
            show_log (bool): If True, the log file will be generated.
            logger (bool): If True, the log file will be generated.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """

        output = directory(self.outputpath) + showbox.removesuffix('.box.pdb')
        command = ['grid', '-i', directory(self.outputpath) + f'{tid}_grid.in', '-o', f'{output}.out', '-t']
        self.perform_subprocess(command, self.outputpath)

        sync_directory(directory(self.outputpath))

    def prepare_gridbox(self):
        """
        This function outlines the steps to prepare a parallel execution of perform_parallel_gridbox method for the DOCK 6 software.
        The analysis is conducted taking into account the number of CPU's specified in the computational architecture.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """
        try:

            self.set_outputpath(f'{self.__base_output_path}/gridbox/')
            showbox_inputpath = f'{self.__base_output_path}/showbox/'

            files = [f for f in os.listdir(directory(showbox_inputpath)) if f.endswith('.box.pdb')]
            args = [(showbox[1], showbox[0]) for showbox in enumerate(files)]

            for tid, showbox in enumerate(files):
                self.generate_docking_script(input_template='src/scripts/dock6/grid.template',
                                            output_script=directory(self.outputpath) + str(tid) + '_grid.in',
                                            receptor_file=directory(self.receptorpath) + self.__pdb_code + '.dockprep.mol2',
                                            box_file=f'{directory(self.__base_output_path)}/showbox/' + showbox,
                                            dock6_path=self.__dock6_path,
                                            score_grid_prefix=directory(self.outputpath) + showbox.removesuffix('.box.pdb'))


            self.process_in_parallel(method_name='perform_parallel_gridbox', args_list=args,process_by_threads=True) if args else None
            if not files:raise ValueError('Nenhuma caixa de busca DOCK6 foi produzida.')
            for showbox in files:
                prefix=Path(self.outputpath)/showbox.removesuffix('.box.pdb')
                for suffix in ('.nrg','.bmp'):
                    file=prefix.with_name(prefix.name+suffix)
                    if not file.is_file() or not file.stat().st_size:raise ValueError('Grid DOCK6 ausente ou vazio: '+file.name)


        except Exception as e:
            self.logger.error(f'during to perform into prepare_gridbox function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            raise

        finally:

            sync_directory(directory(self.outputpath))






    def perform_parallel_minimization(self, ligand:str, tid:int):
        """
        This function outlines the steps to perform the minimization of a ligand using the DOCK 6 software. The analysis is conducted
        in a parallel execution, taking into account the number of CPU's specified in the computational architecture.

        Args:
            ligand (str): The ligand file to perform the minimization.
            tid (int): The identifier of the thread.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """


        output = directory(self.outputpath) + str(tid) + Path(ligand).stem
        command = ['dock6', '-i', directory(self.outputpath) + f'{tid}_min.in', '-o', f'{output}.out']
        self.perform_subprocess(command, self.outputpath)

        sync_directory(directory(self.outputpath))

    def prepare_minimization(self):
        from biomolexplorer.dock6_validation import compound_code
        self.set_outputpath(f'{self.__base_output_path}/energy_min/')
        jobs=[]
        for tid,ligand in enumerate(self._ligand_files(self.ligandpath)):
            code=compound_code(ligand)
            self.generate_docking_script(input_template='src/scripts/dock6/min.template',
                output_script=str(Path(self.outputpath)/f'{tid}_min.in'),
                ligand_atom_file=str(Path(self.ligandpath)/ligand),
                rmsd_reference_filename=str(Path(self.ligandpath)/ligand),
                grid_score_grid_prefix=self._grid_prefix(code),dock6_path=self.__dock6_path,
                ligand_outfile_prefix=str(Path(self.outputpath)/(code+'.lig.min')))
            jobs.append((code,(ligand,tid),Path(self.outputpath)/(code+'.lig.min_scored.mol2')))
        self._run_compound_jobs('perform_parallel_minimization',jobs,'Minimização DOCK6')


    def __identify_residues(self, filename, max_res):
        """
        This function identifies the residues based on the docking analysis. The function reads the output file generated by the DOCK 6 software
        and retrieves the residues with the highest scores. Once the residues are identified, the function returns the residues with the highest
        scores to plot the footprints.

        Args:
            filename (str): The filename to be analyzed.
            max_res (int): The maximum number of residues to be selected.

        Returns:
            resindex_selected (list): The list of residues with the highest scores.
            resindex_remainder (list): The list of residues with the lowest scores.
        """
        try:
            sync_directory(directory(self.outputpath))

            filename = directory(self.outputpath) + filename
            fp_file = open(filename,'r')
            lines = fp_file.readlines()
            fp_file.close()

            num_res = 0
            for line in lines:
                linesplit = line.split()
                if (len(linesplit) == 8):
                    if (linesplit[0] != 'resname'):
                        num_res += 1

            fp_array = [[0 for i in range(2)] for j in range(num_res)]
            for i in range(num_res):
                fp_array[i][0] = i

            count = 0
            for line in lines:
                linesplit = line.split()
                if (len(linesplit) == 8 and linesplit[0] != 'resname'):
                    fp_array[count][1] = max(math.fabs(float(linesplit[2])), math.fabs(float(linesplit[3])), math.fabs(float(linesplit[5])), math.fabs(float(linesplit[6])))
                    count += 1

            fp_array.sort(key=lambda x: x[1])
            resindex_selected = []
            resindex_remainder = []

            max_res = min(max_res, num_res)
            for i in range(max_res):
                resindex_selected.append(fp_array[(num_res-1)-i][0])

            for i in range(num_res - max_res):
                resindex_remainder.append(fp_array[i][0])

            resindex_selected.sort()
            resindex_remainder.sort()
            del fp_array[:][:]

            return resindex_selected, resindex_remainder

        except Exception as e:
            self.logger.error(f'during to perform the __identify_residues function', exc_info=True)
            raise




    def __plot_footprints(self, filename, resindex_selected, resindex_remainder):
        """
        This function plots the footprints based on the docking analysis. The function reads the output file generated by the DOCK 6 software
        and retrieves the residues with the highest scores. Once the residues are identified, the function plots the footprints based on the
        identified residues.

        Args:
            filename (str): The filename to be analyzed.
            resindex_selected (list): The list of residues with the highest scores.
            resindex_remainder (list): The list of residues with the lowest scores.
        """

        sync_directory(directory(self.outputpath))

        data = directory(self.outputpath) + filename + '_footprint_scored.txt'
        footprint = open(data,'r')
        lines = footprint.readlines()
        footprint.close()

        resname = []; resid = []; vdw_ref = []; es_ref = []; vdw_pose = []; es_pose = []
        vdw_score = ""; es_score = ""
        vdw_energy = ""; es_energy = ""
        for line in lines:
            linesplit = line.split()
            if (len(linesplit) == 3):
                if (linesplit[1] in ('vdw_fp:', 'FPS_vdw_fps:')):
                    vdw_score = 'd = '+linesplit[2]
                if (linesplit[1] in ('es_fp:', 'FPS_es_fps:')):
                    es_score = 'd = '+linesplit[2]
                if (linesplit[1] in ('vdw:', 'FPS_vdw_energy:')):
                    vdw_energy = 'vdw = '+linesplit[2]+' kcal/mol'
                if (linesplit[1] in ('es:', 'FPS_es_energy:')):
                    es_energy = 'es = '+linesplit[2]+' kcal/mol'
            if (len(linesplit) == 8):
                if (linesplit[0] != 'resname'):
                    resname.append(linesplit[0])
                    resid.append(linesplit[1])
                    vdw_ref.append(float(linesplit[2]))
                    es_ref.append(float(linesplit[3]))
                    vdw_pose.append(float(linesplit[5]))
                    es_pose.append(float(linesplit[6]))


        resname_selected = []
        vdw_ref_selected = []; es_ref_selected = []; vdw_pose_selected = []; es_pose_selected = []
        for i in (resindex_selected):
            resname_selected.append(resname[i] if resname[i][-1:].isdigit() else resname[i]+resid[i])
            vdw_ref_selected.append(vdw_ref[i])
            es_ref_selected.append(es_ref[i])
            vdw_pose_selected.append(vdw_pose[i])
            es_pose_selected.append(es_pose[i])


        vdw_ref_remainder = 0; es_ref_remainder = 0; vdw_pose_remainder = 0; es_pose_remainder = 0
        for i in (resindex_remainder):
            vdw_ref_remainder += vdw_ref[i]
            es_ref_remainder += es_ref[i]
            vdw_pose_remainder += vdw_pose[i]
            es_pose_remainder += es_pose[i]


        resname_selected.append('REMAIN')
        vdw_ref_selected.append(vdw_ref_remainder)
        es_ref_selected.append(es_ref_remainder)
        vdw_pose_selected.append(vdw_pose_remainder)
        es_pose_selected.append(es_pose_remainder)

        residue = []
        for i in range(len(resname_selected)):
            residue.append(i)

        fig = plt.figure(figsize=(20, 18))
        ax1 = fig.add_subplot(2,1,1)
        ax1.set_title(filename.strip())
        plt.plot(residue, vdw_ref_selected, 'b', linewidth=3)
        plt.plot(residue, vdw_pose_selected, 'r', linewidth=3)
        ax1.set_ylabel('VDW Energy (kcal/mol)')
        ax1.set_xlim(0, len(resname_selected))
        ax1.xaxis.set_major_locator(MultipleLocator(1))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%s'))
        ax1.set_xticks(residue)
        ax1.xaxis.grid(which='major', color='black', linestyle='solid')
        ax1.set_xticklabels(resname_selected, rotation=90)
        ax1.legend(['Reference', 'Pose'])
        ax1.annotate(vdw_score, xy=(0.72,0.12), xycoords='axes fraction', backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        ax1.annotate(vdw_energy, xy=(0.72,0.05), xycoords='axes fraction', backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})

        ax2 = fig.add_subplot(2,1,2)
        plt.plot(residue, es_ref_selected, 'b', linewidth=3)
        plt.plot(residue, es_pose_selected, 'r', linewidth=3)
        ax2.set_ylabel('ES Energy (kcal/mol)')
        ax2.set_xlim(0, len(resname_selected))
        ax2.xaxis.set_major_locator(MultipleLocator(1))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%s'))
        ax2.set_xticks(residue)
        ax2.xaxis.grid(which='major', color='black', linestyle='solid')
        ax2.set_xticklabels(resname_selected, rotation=90)
        ax2.legend(['Reference', 'Pose'])
        ax2.annotate(es_score, xy=(0.72,0.12), xycoords='axes fraction', backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        ax2.annotate(es_energy, xy=(0.72,0.05), xycoords='axes fraction', backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})

        if not os.path.exists(directory(self.outputpath) + 'plots/'):
            os.makedirs(directory(self.outputpath) + 'plots/', exist_ok=True)

        filename = directory(self.outputpath) + 'plots/' + filename + '.pdf'
        plt.savefig(filename)
        plt.close()



    def perform_parallel_footprint(self, receptor:str, ligand:str, tid:int):
        """
        This function outlines the steps to perform the footprint analysis using the DOCK 6 software. The analysis is conducted
        in a parallel execution, taking into account the number of CPU's specified in the computational architecture.

        Args:
            receptor (str): The receptor file to perform the footprint analysis.
            ligand (str): The ligand file to perform the footprint analysis.
            tid (int): The identifier of the thread.
            logger (bool): If True, the log file will be generated.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """

        output = directory(self.outputpath) + receptor + '_' + ligand
        command = ['dock6', '-i', directory(self.outputpath) + f'{tid}_footprint.in', '-o', f'{output}.out']
        self.perform_subprocess(command, self.outputpath)

        sync_directory(directory(self.outputpath))

    def prepare_footprint(self, docked=False):
        from biomolexplorer.dock6_validation import compound_code
        from biomolexplorer.footprints import footprint_receptor
        self.set_outputpath(f'{self.__base_output_path}/footprint/')
        receptor=footprint_receptor(Path(self.receptorpath)/(self.__pdb_code+'.dockprep.mol2'),
                                   Path(self.outputpath)/'footprint_receptor.mol2')
        jobs=[]
        for tid,filename in enumerate(self._ligand_files(self.ligandpath)):
            code=compound_code(filename)
            pose=Path(self.__base_output_path)/(self.__conformer_search_type if docked else 'energy_min')/(code+'_scored.mol2' if docked else code+'.lig.min_scored.mol2')
            self.generate_docking_script(input_template='src/scripts/dock6/footprint.template',
                output_script=str(Path(self.outputpath)/f'{tid}_footprint.in'),ligand_atom_file=str(pose),
                fps_score_footprint_reference_mol2_filename=str(Path(self.ligandpath)/filename),
                fps_score_receptor_filename=str(receptor),dock6_path=self.__dock6_path,
                ligand_outfile_prefix=str(Path(self.outputpath)/code))
            jobs.append((code,(self.__pdb_code,code,tid),Path(self.outputpath)/(code+'_footprint_scored.txt')))
        self._run_compound_jobs('perform_parallel_footprint',jobs,'Footprint DOCK6')


    def plot_footprint_results(self):
        """
        This function plots the footprints based on the docking analysis. The function reads the output file generated by the DOCK 6 software
        and retrieves the residues with the highest scores. Once the residues are identified, the function plots the footprints based on the
        identified residues.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """

        try:

            self.set_outputpath(f'{self.__base_output_path}/footprint/')
            sync_directory(directory(self.outputpath))

            files = [f for f in os.listdir(directory(self.outputpath)) if f.endswith('_footprint_scored.txt')]
            for filename in files:
                from biomolexplorer.footprints import footprint_rows
                code=filename.removesuffix('_footprint_scored.txt');success=False
                with self._compound_guard(code,'Gráfico do footprint DOCK6'):
                    footprint_rows(Path(self.outputpath)/filename)
                    checkpoints=self._checkpoints() if hasattr(self,'runtime_path') else None
                    plot=Path(self.outputpath)/'plots'/(code+'.pdf')
                    if checkpoints:
                        from biomolexplorer.calculation_cache import identity
                        from biomolexplorer.stage_cache import file_digest
                        key=identity(dict(footprint=file_digest(Path(self.outputpath)/filename),
                                          max_residues=self.__max_residues,adapter=checkpoints.adapter,kind='plot'))
                        if checkpoints.cache.restore(key,{'plot':plot}) is not None:continue
                    resindex_selected, resindex_remainder = self.__identify_residues(filename, self.__max_residues)
                    self.__plot_footprints(code, resindex_selected, resindex_remainder)
                    if checkpoints:checkpoints.cache.remember(key,{'plot':plot},dict(compound=code,step='footprint-plot'))
                    success=True
                if not success:self.excluded_compounds.add(code)
            self.failure_report.write(self.export_path)

        except Exception as e:
            self.logger.error(f'during to perform the plot_footprint_results function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            raise

        finally:

            sync_directory(directory(self.outputpath))




    def perform_parallel_docking(self, tid:int):
        """
        This function outlines the steps to perform the docking analysis using the DOCK 6 software. The analysis is conducted
        in a parallel execution, taking into account the number of CPU's specified in the computational architecture.

        Args:
            tid (int): The identifier of the thread.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """

        command = ['dock6', '-i', directory(self.outputpath) + f'{tid}_docking.in',
                   '-o', directory(self.outputpath) + f'{tid}_docking.out']
        self.perform_subprocess(command, self.outputpath)

        sync_directory(directory(self.outputpath))

    def perform_dock6_evaluation(self):
        from biomolexplorer.dock6_validation import compound_code
        self.set_outputpath(f'{self.__base_output_path}/{self.__conformer_search_type}/')
        jobs=[]
        for tid,ligand in enumerate(self._ligand_files(Path(self.__base_output_path)/'energy_min')):
            code=compound_code(ligand)
            spheres='binding_site' if getattr(self,'binding_site_center',None) is not None else code
            self.generate_docking_script(input_template='src/scripts/dock6/docking.template',
                output_script=str(Path(self.outputpath)/f'{tid}_docking.in'),
                conformer_search_type=self.__conformer_search_type,
                ligand_atom_file=str(Path(self.__base_output_path)/'energy_min'/ligand),
                rmsd_reference_filename=str(Path(self.__base_output_path)/'energy_min'/ligand),
                receptor_site_file=str(Path(self.__base_output_path)/'surface/Molecules'/(spheres+'.sph')),
                ligand_sphere_file=str(Path(self.__base_output_path)/'surface/Molecules'/(spheres+'.sph')),
                grid_score_grid_prefix=self._grid_prefix(code),dock6_path=self.__dock6_path,
                ligand_outfile_prefix=str(Path(self.outputpath)/code))
            jobs.append((code,(tid,),Path(self.outputpath)/(code+'_scored.mol2')))
        self._run_compound_jobs('perform_parallel_docking',jobs,'Docking DOCK6')
