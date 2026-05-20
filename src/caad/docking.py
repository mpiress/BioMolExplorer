"""!

PhD research 2023~2026

@title 
    Drug Hunter: a general framework based on a target-ligand strategy to investigate the 
    physicochemical potential of compounds

@info
    A general crawler to look at targets and bioactivities molecules on the ChEMBL database.

@authors 
   - Michel Pires da Silva (michel@cefetmg.br / Academic student)
   - Alisson Marques da Silva (alisson@cefetmg.br / Co Advisor)
   - Alex Gutterres Taranto   (taranto@ufsj.edu.br / Advisor)

@date 2023-2026

@copyright GNU Public License

@cond GNU_PUBLIC_LICENSE
    Drug Hunter is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    Drug Hunter is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
@endcond

"""
#----------------------------------------------------------------------------------------------
#Configure PYTHONPATH to perform execution using the project classes
import sys 
sys.path.append("src")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
import warnings
# Desabilita todos os avisos
warnings.filterwarnings("ignore")

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
        self.nprocess       = os.cpu_count() - 2
        self.centers        = None
        self.logpath        = self.path + '/logs/'
        self.logger         = LoggerManager.get_logger(self.__class__.__name__, log_file='logs/docking.log')
        self.mol_filename   = mol_filename
        

        self.set_outputpath(output_path) if output_path != None else None
        
        


    def set_ligandpath(self, path) -> None:
       if not os.path.exists(self.path + path):
           print(f'[ERROR]: The ligand path {path} does not exist!')
           exit(1)

       self.ligandpath = path
    
    
    def set_receptorpath(self, path) -> None:
       if not os.path.exists(self.path + path):
           print(f'[ERROR]: The receptor path {path} does not exist!')
           exit(1)
       
       self.receptorpath = path
       
       
    def set_complexpath(self, path) -> None:
       if not os.path.exists(self.path + path):
           print(f'[ERROR]: The complex path {path} does not exist!')
           exit(1)
       
       self.complexpath = path
       
       
    def set_outputpath(self, path) -> None:
        if not os.path.exists(self.path + path):
            os.makedirs(self.path + path, exist_ok=True)
        
        if not os.path.exists(self.path + path):
            print(f'[ERROR]: The output path {path} can not be created!')
            exit(1)

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
            with open(input_template, 'r') as template_file:
                template_content = template_file.read()

            config_content = template_content.format(**kwargs)

            with open(output_script, 'w') as output_file:
                output_file.write(config_content)

        except Exception as e:
            self.logger.error(f'during to perform {input_template} -> {output_script} in generate_docking_script function', exc_info=True)
        
        finally:
            time.sleep(1)
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
           



    def perform_subprocess(self, command:str, local_path=None, shell=True, check=True) -> bool:
        """
        This function performs a subprocess execution of a given command. 
        """
        
        try:

            if local_path != None:
                tmp = subprocess.run(
                    command, cwd=local_path[1:], 
                    shell=shell, check=check, start_new_session=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                
            else:
                tmp = subprocess.run(
                    command, 
                    shell=shell, check=check, start_new_session=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )


            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f'STDERR: {e.stdout}', exc_info=True)
            self.logger.error(f'STDERR: {e.stderr}', exc_info=True)
            return False
        
        finally:
            if local_path != None:
                dir_fd = os.open(local_path[1:], os.O_DIRECTORY)
                os.fsync(dir_fd)
            else:
                dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
                os.fsync(dir_fd)
            time.sleep(1)
        
        
        

    
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
            command = f'chimera --nogui --silent {filename}'
            self.perform_subprocess(command, self.outputpath)
 
        except Exception as e:
            self.logger.error(f'during to perform {filename} in prepare_on_chimera function', exc_info=True)
        
        finally:
            time.sleep(1)
            os.remove(self.outputpath[1:] + filename) if os.path.isfile(self.outputpath[1:] + filename) else None
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
                



    
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
            
            input_file = self.path + self.outputpath + inputfile 
            
            if input_format != 'smi':
                command = f'obabel -i {input_format} "{input_file}" -o {output_format} -O "{outputfile}" '
            else:
                command = f'obabel -:"{inputfile}" -o {output_format} -O "{outputfile}" --gen3D '
            
            
            for param, value in params: 
                command += f'-{param} {value} ' if value else f'-{param} '
            
            self.perform_subprocess(command, self.outputpath)
   
        except Exception as e:
            self.logger.error(f'during to perform {inputfile} to {outputfile} converter in prepare_on_obabel function', exc_info=True)

        finally:
            time.sleep(1)
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            
            
            


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

            file = self.path + self.outputpath + inputfile
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
            return self.centers[center_idx].values.tolist() if center_idx in self.centers else None
            
        except Exception as e:
            self.logger.error(f'during to perform {receptor} -> {ligand }in retrieve_centerofmass_dataset function', exc_info=True)



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
               
            
    
    def prepare_for_docking(self, pdb_codes:list, charge_type:str, pH:float, redefine_centerofmass:bool) -> bool:
        """
        This function prepares the input files for docking using Chimera and Open Babel software. The input files are prepared in a specific
        format required to perform docking with different docking software. Templates for the input files are available in the src/scripts folder.
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
            
            centers   = {} 
            complexes = []
            receptors = []
            ligands   = []


            for pdb in pdb_codes:
                receptor = pdb[0]
                ligand   = pdb[1]
                resnum   = pdb[2] 
                chain    = pdb[3]
                chain_id = ',.'.join(chain) if len(chain) > 1 else chain
                
                if f'prepare_complex_{receptor}_{chain}.com' not in complexes:
                    self.generate_docking_script(input_template='src/scripts/chimera/prepare_complex.template',
                                                output_script=self.outputpath[1:] + f'prepare_complex_{receptor}_{chain}.com',
                                                pdb_code=self.path + self.complexpath + receptor,
                                                input_complex=self.path + self.complexpath + receptor,
                                                chain=chain_id,
                                                output_complex=f'{receptor}_{chain}')
                    complexes.append(f'prepare_complex_{receptor}_{chain}.com')
                
                
                if f'prepare_receptor_{receptor}_{chain}.com' not in receptors:
                    self.generate_docking_script(input_template='src/scripts/chimera/prepare_receptor.template',
                                                output_script=self.outputpath[1:] + f'prepare_receptor_{receptor}_{chain}.com',
                                                input_complex=self.path + self.outputpath + f'{receptor}_{chain}' + '.complex',
                                                receptor=self.path + self.outputpath + f'{receptor}_{chain}')
                    receptors.append(f'prepare_receptor_{receptor}_{chain}.com')
            
                
                extention = f'{receptor}_{ligand}_{resnum}{chain}'
                self.generate_docking_script(input_template='src/scripts/chimera/prepare_ligand.template',
                                             output_script=self.outputpath[1:] + f'prepare_ligand_{extention}.com',
                                             input_complex=self.path + self.outputpath + f'{receptor}_{chain}' + '.complex',
                                             resnum=resnum,
                                             chain=chain[0],
                                             charge_type=charge_type,
                                             input_ligand=self.path + self.outputpath + f'{extention}',
                                             output_ligand=self.path + self.outputpath + f'{extention}')
                ligands.append(f'prepare_ligand_{extention}.com')
               
            

            args = [(file,) for file in complexes]
            self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)

            args = [(file,) for file in receptors]
            self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)
            
            for pdb in pdb_codes:
                self.prepare_on_obabel(f'{pdb[0]}_{chain}.dockprep.mol2', f'{pdb[0]}_{chain}.dockprep.pdbqt', [('p',str(pH)), ('xr','')], input_format='mol2', output_format='pdbqt')
            
            args = [(file,) for file in ligands]
            self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)  
            
            tmp_codes = []
            for pdb in pdb_codes:
                key = f'{pdb[0]}_{pdb[1]}_{pdb[2]}{chain}'
                self.prepare_on_obabel(f'{key}.lig.mol2', f'{key}.lig.pdbqt', [('p',str(pH))], input_format='mol2', output_format='pdbqt')
                if redefine_centerofmass:
                    center = self.calculate_ligand_centerofmass(f'{key}.lig.pdb', pdb[1])
                    if center != None:
                        centers[key] = center
                        tmp_codes.append(tuple(pdb)) 

            
            if redefine_centerofmass:
                f1 = fileHandling(input_path=self.outputpath, output_path=self.outputpath)
                tmp = f1.csv_to_dataframe('centers')
                tmp = tmp.to_dict(orient='list')
                tmp.update(centers)
                df = DataFrame(tmp)
                f1.dataframe_to_csv('centers', df)

            return tmp_codes
            

        except Exception as e:
             return None
        
        finally:
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(1)
                


    def recover_better_conforms_of_vina(self, charge_type:str, filename:Optional[list]=None, molecules_dataset:Optional[str]=None):  
        """
        This function recover better conformations for each ligand performed by AutoDock Vina software. The input files are prepared in a specific
        format required to perform docking with different docking software. Templates for the input files are available in the src/scripts folder.
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

            files = [f for f in os.listdir(self.ligandpath[1:]) if f.endswith('.pdbqt')] if filename == None else filename
            conv      = MolConverter(input_path=self.ligandpath, output_path=self.outputpath)
            explorer  = MolExplorer(input_path=self.ligandpath, output_path=self.outputpath)
            
            for input_file in files:

                pdb_str_content = conv.extract_pdb_to_pdbqt(input_file, start_index='MODEL 1', end_index='MODEL 2', pdb_filename=input_file.rsplit('.')[0]+'.pdb')
               
                
                if not explorer.is_fragmented(pdb_str_content):
                    
                    self.generate_docking_script(input_template='src/scripts/chimera/prepare_better_conform.template',
                                             output_script=self.outputpath[1:] + f'prepare_better_conform_{input_file}.com',
                                             ligand=input_file.rsplit(".")[0], charge_type=charge_type)

            
            args = [(f'prepare_better_conform_{input_file}.com',) for input_file in files]
            self.process_in_parallel(method_name='prepare_on_chimera', args_list=args)
            

            files = [f for f in os.listdir(self.outputpath[1:]) if f.endswith('.lig.mol2')]
            for file in files:
                with open(self.outputpath[1:] + file, 'r') as fp:
                    lines = ''.join(fp.readlines())
                    if lines.find('nan') >= 0:
                        os.remove(self.outputpath[1:] + file)


            if molecules_dataset:
                fx = fileHandling(input_path=molecules_dataset, output_path=molecules_dataset)
                files = [f.replace('.lig.mol2','').split('_')[1] for f in os.listdir(self.outputpath[1:]) if f.endswith('.lig.mol2')]
                df = fx.csv_to_dataframe(self.mol_filename)
                data = df[df['molecule_chembl_id'].isin(files)]
                fx.dataframe_to_csv(self.mol_filename, data)
                data = set(files) - set(data['molecule_chembl_id'].tolist())
                path = self.path + self.outputpath
                files = [os.remove(path + f) for f in os.listdir(path) if f.replace('.lig.mol2','').split('_')[1] in data]
                

        except Exception as e:
            self.logger.error(f'during to perform {filename} in recover_better_conformation function', exc_info=True)

        finally:
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(1)
            
            pdbs = [f for f in os.listdir(self.outputpath[1:]) if f.endswith('.pdb')]
            [os.remove(self.outputpath[1:] + file) for file in pdbs]
            
            mols = [f for f in os.listdir(self.outputpath[1:]) if f.endswith('(2).lig.mol2')]
            [os.remove(self.outputpath[1] + file) for file in mols]
            
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(1)
            
        

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
        try:

            f1 = fileHandling(input_path=self.ligandpath, output_path=self.ligandpath)
            df = f1.csv_to_dataframe(self.mol_filename)
            df['molecule_chembl_id'] = df['molecule_chembl_id'].astype(str)
            data = df[['molecule_chembl_id', 'canonical_smiles']].to_records(index=False)
            
            for chemblid, smiles in data:
                self.prepare_on_obabel(inputfile=smiles, outputfile=chemblid + '.lig.pdbqt', input_format="smi", output_format='pdbqt', params=[('p',str(pH))])
    

        except Exception as e:
            self.logger.error(f'during to perform {chemblid} -> {smiles} in prepare_compounds function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)

        finally:
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(1)

            mols = [f.split('.lig.pdbqt')[0] for f in os.listdir(self.outputpath[1:]) if f.endswith('.lig.pdbqt')]
            df = df[df['molecule_chembl_id'].isin(mols)]
            f1.dataframe_to_csv(self.mol_filename, df)

            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(1)



    def perform_vina_evaluation(self):
        """
        This function performs the AutoDock Vina evaluation for each prepared ligand in the ligand path folder. The assessment
        is performed sequentially. The results are stored in the output path folder in PDBQT files named according to the
        identifications reported in molecules.csv.

        Raises:
            Exception: If any error occurs during the evaluation of the prepared ligands for docking, a log file will
            be created and posted in the log folder.

        """
        
        files_to_perform = [f for f in os.listdir(self.outputpath[1:]) if f.endswith('.vina')] 
            
        for file in files_to_perform:
            command = f'vina --config {file}'
            self.perform_subprocess(command, self.outputpath)
            
       
    
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
            idx_to_remove = []

            for idx, (receptor, ligand, resnum, chain) in enumerate(pdb_codes):
                composite = f'{receptor}_{ligand}_{resnum}{chain}'
                
                center = self.retrieve_centerofmass_dataset(self.ligandpath, receptor, ligand, resnum, chain)
                if center == None:
                    idx_to_remove.append(idx)
                    continue

                self.generate_docking_script(input_template='src/scripts/vina/config.template',
                                            output_script=self.outputpath[1:] + f'{composite}.vina',
                                            receptor=self.path + self.receptorpath + f'{receptor}_{chain}.dockprep.pdbqt',
                                            ligand=self.path + self.ligandpath + f'{composite}' + '.lig.pdbqt',
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
            
            pdb_codes = np.delete(pdb_codes, idx_to_remove)
            for receptor, ligand, resnum, chain in pdb_codes:
                composite = f'{receptor}_{ligand}_{resnum}{chain}'
                iligand  = self.path + self.ligandpath + f'{composite}' + '.lig.pdbqt'
                vina_model = self.path + self.outputpath + f'{composite}' + '.lig.pdbqt'
                if os.path.isfile(iligand) and os.path.isfile(vina_model):
                    results.append((f'{receptor}', f'{ligand}', f'{resnum}', f'{chain}', desc.calcRMSD(iligand, vina_model)))

            
            rmsd = DataFrame(results, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RMSD'])
            pdb_codes = DataFrame(self.__pdb_codes, columns=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN', 'RESOLUTION'])
            pdb_codes['RESNUM'] = pdb_codes['RESNUM'].astype(str)
            pdb_codes = pdb_codes.merge(rmsd, on=['PDB_CODE', 'LIGAND', 'RESNUM', 'CHAIN'], how='left')
            
            f1   = fileHandling(output_path=self.complexpath)
            f1.dataframe_to_csv('pdb_codes', pdb_codes)
                    
               
        except Exception as e:
            self.logger.error('during to perform the docking function', exc_info=True)



    def docking(self, base_selected_mols:str):

        """
        This function performs the docking analysis using the AutoDock Vina software. The analysis is conducted sequentially, considering the number of processes
        available in the molecules.csv. The search box size is set to [24,24,24], the exhaustiveness is set to 20, and the number of result modes is set to 10. The parameters
        can be modified based on user conditions, and results are stored in a series of PDBQT files, named according to the identifiers reported in molecules.csv.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.

        """

        try:
            
            molecules = [f.rsplit('.lig.pdbqt')[0] for f in os.listdir(self.ligandpath[1:]) if f.endswith('.lig.pdbqt')]
            
            self.__pdb_codes   = [(pdb[0], pdb[1], pdb[2], pdb[3]) for pdb in self.__pdb_codes]
            
            for receptor, ligand, resnum, chain in self.__pdb_codes:
                center    = self.retrieve_centerofmass_dataset(self.__centerofmasspath, receptor, ligand, resnum, chain)
                tmp       = [f.replace('.lig.pdbqt','').replace(f'{receptor}_', '') for f in os.listdir(self.outputpath[1:]) if f.endswith('.lig.pdbqt')]
                molecules = set(molecules) - set(tmp)
               
                for mol in molecules:
                     self.generate_docking_script(input_template='src/scripts/vina/config.template',
                                            output_script=self.outputpath[1:] + f'{receptor}_{mol}.vina',
                                            receptor=self.path + self.receptorpath + f'{receptor}_{chain}' + '.dockprep.pdbqt',
                                            ligand=self.path + self.ligandpath + mol + '.lig.pdbqt',
                                            center_x=center[0],
                                            center_y=center[1],
                                            center_z=center[2],
                                            size_x=self.__sizeof_box[0],
                                            size_y=self.__sizeof_box[1],
                                            size_z=self.__sizeof_box[2],
                                            out=f'{receptor}_{mol}.lig.pdbqt',
                                            exhaustiveness=self.__exhaustiveness,
                                            num_modes=self.__num_modes)
                     
                     time.sleep(1)
                     dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
                     os.fsync(dir_fd)
                     time.sleep(1)
                     
                     command  = f'vina --config {receptor}_{mol}.vina'
                     validate = self.perform_subprocess(command, self.outputpath)
                     
                     if not validate:
                         f1 = fileHandling(input_path=base_selected_mols, output_path=base_selected_mols)
                         df = f1.csv_to_dataframe(self.mol_filename)
                         df = df[df['molecule_chembl_id'] != mol]
                         f1.dataframe_to_csv(self.mol_filename, df)
                         os.remove(self.path + self.ligandpath + mol + '.lig.pdbqt') if os.path.isfile(self.ligandpath[1:] + mol + '.lig.pdbqt') else None
                         

        except Exception as e:
            self.logger.error('during to perform the docking function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)

        finally:
            time.sleep(1)
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(1)

            [os.remove(self.outputpath[1:] + file) for file in os.listdir(self.outputpath[1:]) if file.endswith('.vina')]
            [os.remove(self.outputpath[1:] + file) for file in os.listdir(self.outputpath[1:]) if file.endswith('(2).pdbqt')]
            self.centers = None
                    



class Dock6(Docking):
    
    
    def __init__(self, dock6_path:Optional[str]='', ligand_input_path:Optional[str]=None, receptor_input_path:Optional[str]=None,
                 base_output_path:Optional[str]=None, pdb_code:Optional[str]=None, density:Optional[float]=0.5, 
                 radius:Optional[float]=1.4, distance:Optional[float]=10.0, max_residues:Optional[int]=50,
                 conformer_search_type:Optional[Literal['flex', 'rigid']] = 'flex', mol_filename:Optional[str]='molecules',) -> None:
        
        super().__init__(ligand_input_path=ligand_input_path,
                         receptor_input_path=receptor_input_path,
                         output_path=f'{base_output_path}/',
                         mol_filename=mol_filename)
        
        self.__dock6_path            = dock6_path
        self.__base_output_path      = base_output_path  
        self.__pdb_code              = pdb_code
        self.__density               = density
        self.__radius                = radius
        self.__distance              = distance
        self.__max_residues          = max_residues
        self.__conformer_search_type = conformer_search_type


     
    def prepare_surface(self) -> None:
        
        """
        The following function outlines the steps to prepare the surface of a pdb_code using the DOCK 6 software. 
        As default, the density is set to 0.5 and the radius to 1.4, but, if necessary, these parameters can be modified.

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.

        """

        try:

            self.set_outputpath(f'{self.__base_output_path}/surface/' )
            
            os.remove(self.outputpath[1:] + self.__pdb_code +".dms") if os.path.exists(self.outputpath[1:] + self.__pdb_code +".dms") else None
            os.remove(self.outputpath[1:] + self.__pdb_code +".sph") if os.path.exists(self.outputpath[1:] + self.__pdb_code +".sph") else None
            
            input = self.receptorpath[1:] + self.__pdb_code
            output = self.outputpath[1:] + self.__pdb_code

            command = f'dms "{input}.noH.pdb" -d {self.__density} -n -w {self.__radius} -v -o {output}.dms'
            self.perform_subprocess(command)
            
            self.generate_docking_script(input_template='src/scripts/dock6/INSPH.template',
                                         output_script=self.outputpath[1:]+'INSPH',
                                         receptor=self.__pdb_code)
            
            command = f'sphgen -i INSPH -o OUTSPH'
            self.perform_subprocess(command, self.outputpath)
            
            #prepare the surface selectors
            self.set_outputpath(f'{self.__base_output_path}/surface/Molecules/')
            selector = self.path + f'{self.__base_output_path}/surface/' + self.__pdb_code + '.sph'
            
            files = [f for f in os.listdir(self.ligandpath[1:]) if f.endswith('.mol2')]
            for ligand in files:
                command = f'sphere_selector {selector} {self.path + self.ligandpath + ligand} {self.__distance}'
                self.perform_subprocess(command, self.outputpath)
                os.rename(self.outputpath[1:] + 'selected_spheres.sph', self.outputpath[1:] + ligand.rsplit('.')[0] + '.sph')
    
        except Exception as e:
            self.logger.error(f'during to perform into prepare_surface function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            
        finally:
            time.sleep(1)
            output = f'{self.__base_output_path[1:]}/surface/' 

            os.remove(output + "temp1.ms") if os.path.exists(output + "temp1.ms") else None
            os.remove(output + "temp2.sph") if os.path.exists(output + "temp2.sph") else None
            os.remove(output + "temp3.atc") if os.path.exists(output + "temp3.atc") else None
            os.remove(output + "OUTSPH") if os.path.exists(output + "OUTSPH") else None
            os.remove(output + 'INSPH') if os.path.exists(output + 'INSPH') else None

            dir_fd = os.open(output, os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(1)
            


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
            
            files = [f for f in os.listdir(surface_inputpath[1:]) if f.endswith('.sph')]
            for ligand in files:
                
                self.generate_docking_script(input_template='src/scripts/dock6/showbox.template',
                                             output_script=self.outputpath[1:] + f'{ligand.split(".")[0]}.in',
                                             in_surface='../surface/Molecules/' + ligand,
                                             out_surface=ligand.rsplit('.')[0]+'.box.pdb')
    
            
            for ligand in files:
                command = f'showbox < {ligand.split(".")[0]}.in'
                self.perform_subprocess(command, self.outputpath)

            
        except Exception as e:
            self.logger.error(f'during to perform the prepare_showbox function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
        
        finally:
            time.sleep(1)
            files = [f for f in os.listdir(self.outputpath[1:]) if f.endswith('.in')]
            [os.remove(self.outputpath[1:] + ligand) for ligand in files if os.path.exists(self.outputpath[1:] + ligand)]
            time.sleep(1)
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd) 
            

            
            
            
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
        
        output = self.outputpath[1:] + showbox.rsplit('.')[0]
        command = f'grid -i {self.outputpath[1:] + str(tid)}_grid.in -o {output}.out -t'
        self.perform_subprocess(command)

        dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
        os.fsync(dir_fd)
        time.sleep(1)
               


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
            
            files = [f for f in os.listdir(showbox_inputpath[1:]) if f.endswith('.box.pdb')]
            args = [(showbox[1], showbox[0]) for showbox in enumerate(files)]

            for tid, showbox in enumerate(files):
                self.generate_docking_script(input_template='src/scripts/dock6/grid.template',
                                            output_script=self.outputpath[1:] + str(tid) + '_grid.in',
                                            receptor_file=self.receptorpath[1:] + self.__pdb_code + '.dockprep.mol2',
                                            box_file=f'{self.__base_output_path[1:]}/showbox/' + showbox,
                                            dock6_path=self.__dock6_path,
                                            score_grid_prefix=self.outputpath[1:] + showbox.rsplit('.')[0])


            self.process_in_parallel(method_name='perform_parallel_gridbox', args_list=args) if args else None


        except Exception as e:
            self.logger.error(f'during to perform into prepare_gridbox function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
            
        finally:
            time.sleep(1)
            [os.remove(self.outputpath[1:] + f) for f in os.listdir(self.outputpath[1:]) if f.endswith('_grid.in')]
            [os.remove(self.outputpath[1:] + f) for f in os.listdir(self.outputpath[1:]) if f.endswith('.out')] 
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            
            
            



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
       

        output = self.outputpath[1:] + str(tid) + ligand.rsplit('.')[0]
        command = f'dock6 -i {self.outputpath[1:] + str(tid)}_min.in -o {output}.out'
        self.perform_subprocess(command)

        dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
        os.fsync(dir_fd)
        time.sleep(1)
                  
        
    
    
    def prepare_minimization(self):
        """
        This function outlines the steps to prepare a parallel execution of perform_parallel_minimization method for the DOCK 6 software. 
        The analysis is conducted taking into account the number of CPU's specified in the computational architecture. 

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """ 
        
        try:
            self.set_outputpath(f'{self.__base_output_path}/energy_min/') 

            files = [f for f in os.listdir(f'{self.path+self.ligandpath}') if f.endswith('.mol2')]
            args = [(ligand[1], ligand[0]) for ligand in enumerate(files)]

            for tid, ligand in enumerate(files):
                self.generate_docking_script(input_template='src/scripts/dock6/min.template',
                                        output_script=self.outputpath[1:] + str(tid) +"_min.in",
                                        ligand_atom_file=f'{self.ligandpath[1:]}' + ligand,
                                        rmsd_reference_filename=f'{self.ligandpath[1:]}' + ligand,
                                        grid_score_grid_prefix=f'{self.__base_output_path[1:]}/gridbox/' + ligand.rsplit('.')[0],
                                        dock6_path=self.__dock6_path,
                                        ligand_outfile_prefix=self.outputpath[1:] + ligand.rsplit('.')[0]+'.lig.min')

            
            
            self.process_in_parallel(method_name='perform_parallel_minimization', args_list=args) if args else None    
            

        except Exception as e:
            self.logger.error(f'during to perform into prepare_minimization function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)

        finally:
            time.sleep(5)
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            if os.path.isdir(self.path + self.outputpath):
                [os.remove(self.outputpath[1:] + f) for f in os.listdir(self.outputpath[1:]) if f.endswith('_min.in')]
                [os.remove(self.outputpath[1:] + f) for f in os.listdir(self.outputpath[1:]) if f.endswith('.out')] 
                dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
                os.fsync(dir_fd)
            
             



    
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
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)

            filename = self.path + self.outputpath + filename
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
        
        dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
        os.fsync(dir_fd)

        data = self.path + self.outputpath + filename + '_footprint_scored.txt'
        footprint = open(data,'r')
        lines = footprint.readlines()
        footprint.close()

        resname = []; resid = []; vdw_ref = []; es_ref = []; vdw_pose = []; es_pose = []
        vdw_score = ""; es_score = ""
        vdw_energy = ""; es_energy = ""
        for line in lines:
            linesplit = line.split()
            if (len(linesplit) == 3):
                if (linesplit[1] == 'vdw_fp:'):
                    vdw_score = 'd = '+linesplit[2]
                if (linesplit[1] ==  'es_fp:'):
                    es_score = 'd = '+linesplit[2]
                if (linesplit[1] == 'vdw:'):
                    vdw_energy = 'vdw = '+linesplit[2]+' kcal/mol'
                if (linesplit[1] == 'es:'):
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
            resname_selected.append(resname[i]+resid[i])
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
        ax1.set_ylabel('VDW Energy')
        ax1.set_ylim(-10, 5)
        ax1.set_xlim(0, len(resname_selected))
        ax1.xaxis.set_major_locator(MultipleLocator(1))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%s'))
        ax1.set_xticks(residue)
        ax1.xaxis.grid(which='major', color='black', linestyle='solid')
        ax1.set_xticklabels(resname_selected, rotation=90)
        ax1.legend(['Reference', 'Pose'])
        ax1.annotate(vdw_score, xy=(37,-8), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        ax1.annotate(vdw_energy, xy=(37,-9), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        
        ax2 = fig.add_subplot(2,1,2)
        plt.plot(residue, es_ref_selected, 'b', linewidth=3)
        plt.plot(residue, es_pose_selected, 'r', linewidth=3)
        ax2.set_ylabel('ES Energy')
        ax2.set_ylim(-10, 5)
        ax2.set_xlim(0, len(resname_selected))
        ax2.xaxis.set_major_locator(MultipleLocator(1))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%s'))
        ax2.set_xticks(residue)
        ax2.xaxis.grid(which='major', color='black', linestyle='solid')
        ax2.set_xticklabels(resname_selected, rotation=90)
        ax2.legend(['Reference', 'Pose'])
        ax2.annotate(es_score, xy=(37,-8), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        ax2.annotate(es_energy, xy=(37,-9), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        
        if not os.path.exists(self.outputpath[1:] + 'plots/'):
            os.makedirs(self.outputpath[1:] + 'plots/', exist_ok=True)
            
        filename = self.outputpath[1:] + 'plots/' + filename + '.pdf'
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
            
        output = self.outputpath[1:] + receptor + '_' + ligand
        command = f'dock6 -i {self.outputpath[1:] + str(tid)}_footprint.in -o {output}.out'
        self.perform_subprocess(command)
            
        dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
        os.fsync(dir_fd)
        time.sleep(1)
               
        
    
    def prepare_footprint(self):
        """
        This function outlines the steps to prepare a parallel execution of perform_parallel_footprint method for the DOCK 6 software. 
        The analysis is conducted taking into account the number of CPU's specified in the computational architecture. 

        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """

        try:
            self.set_outputpath(f'{self.__base_output_path}/footprint/')

            files = [f.split('.')[0] for f in os.listdir(f'{self.ligandpath[1:]}') if f.endswith('.mol2')] 
            
            args = [(self.__pdb_code, ligand[1], ligand[0]) for ligand in enumerate(files)]
            
            for tid, ligand in enumerate(files):
                self.generate_docking_script(input_template='src/scripts/dock6/footprint.template',
                                        output_script=self.outputpath[1:] + str(tid) +"_footprint.in",
                                        ligand_atom_file=f'{self.__base_output_path[1:]}/energy_min/' + ligand + '.lig.min_scored.mol2',
                                        fps_score_footprint_reference_mol2_filename=f'{self.ligandpath[1:]}' + ligand + '.lig.mol2',
                                        fps_score_receptor_filename=self.receptorpath[1:] + self.__pdb_code + '.dockprep.mol2',
                                        dock6_path=self.__dock6_path,
                                        ligand_outfile_prefix=self.outputpath[1:] + ligand)


            self.process_in_parallel(method_name='perform_parallel_footprint', args_list=args) if args else None
        

        except Exception as e:
            self.logger.error(f'during to perform into prepare_footprint function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)


        finally:
            time.sleep(1)
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            [os.remove(self.outputpath[1:] + f) for f in os.listdir(f'{self.outputpath[1:]}') if f.endswith('.in')]
            [os.remove(self.outputpath[1:] + f) for f in os.listdir(f'{self.outputpath[1:]}') if f.endswith('.out')]
            


    
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
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            time.sleep(5)

            files = [f for f in os.listdir(self.outputpath[1:]) if f.endswith('_footprint_scored.txt')]
            for filename in files:
                resindex_selected, resindex_remainder = self.__identify_residues(filename, self.__max_residues)
                self.__plot_footprints(filename.replace('_footprint_scored.txt', ''), resindex_selected, resindex_remainder)
        
        except Exception as e:
            self.logger.error(f'during to perform the plot_footprint_results function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
        
        finally:
            time.sleep(1)
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
            
            
    
       
    def perform_parallel_docking(self, tid:int):
        """
        This function outlines the steps to perform the docking analysis using the DOCK 6 software. The analysis is conducted 
        in a parallel execution, taking into account the number of CPU's specified in the computational architecture. 

        Args:
            tid (int): The identifier of the thread.
            
        Raises:
            Exception: If any error occurs during the docking analysis, a log file will be created and posted in the log folder.
        """
        
        command = f'dock6 -i {self.outputpath[1:]+str(tid)}_docking.in'
        self.perform_subprocess(command)

        dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
        os.fsync(dir_fd)
        time.sleep(1)

        
    
        
    def perform_dock6_evaluation(self):
        
        try:
            self.set_outputpath(f'{self.__base_output_path}/flex/')
            
            files = [f for f in os.listdir(f'{self.__base_output_path[1:]}/energy_min/') if f.endswith('.mol2')]
            args = [(ligand[0],) for ligand in enumerate(files)]

            for tid, ligand in enumerate(files):
                self.generate_docking_script(input_template='src/scripts/dock6/docking.template',
                                        output_script=self.outputpath[1:] + str(tid) + '_docking.in',
                                        conformer_search_type=self.__conformer_search_type,
                                        ligand_atom_file=f'{self.__base_output_path[1:]}/energy_min/' + ligand,
                                        rmsd_reference_filename=f'{self.__base_output_path[1:]}/energy_min/' + ligand,
                                        receptor_site_file=f'{self.__base_output_path[1:]}/surface/Molecules/{ligand.split('.')[0]}.sph', 
                                        ligand_sphere_file=f'{self.__base_output_path[1:]}/surface/Molecules/{ligand.split('.')[0]}.sph',
                                        grid_score_grid_prefix=f'{self.__base_output_path[1:]}/gridbox/' + ligand.split('.')[0],
                                        dock6_path=self.__dock6_path,
                                        ligand_outfile_prefix=self.outputpath[1:] + ligand.split('.')[0])
                
            
            self.process_in_parallel(method_name='perform_parallel_docking', args_list=args) if args else None
        
        except Exception as e:
            self.logger.error(f'during to perform the perform_dock6_evaluation function', exc_info=True)
            self.logger.error(f'STDERR: {e}', exc_info=True)
        
        finally:
            time.sleep(1)
            [os.remove(self.outputpath[1:] + f) for f in os.listdir(self.outputpath[1:]) if f.endswith('.in')]
            dir_fd = os.open(self.outputpath[1:], os.O_DIRECTORY)
            os.fsync(dir_fd)
           
        
        
        