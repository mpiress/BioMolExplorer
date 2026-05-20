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
from wrappers.docking import perform_consensus
#----------------------------------------------------------------------------------------------


if __name__ == "__main__":

    #----------------------------------------------------------------------------------------------
    # APÓS IDENTIFICAR O MELHOR PAR RECEPTOR LIGANTE É PRECISO VERIFICAR SE HÁ LOOPS FALTANTES
    # EM CASO AFIRMATIVO, PARE AQUI A EXECUÇÃO, FAÇA AS CORREÇÕES E DEPOIS EXECUTE AS CHAMADAS 
    # A BAIXO NOS PARES IDENTIFICADOS PARA REFATORAR SUAS REPRESENTAÇÕES NO DATASET.
    #----------------------------------------------------------------------------------------------
    
    #perform_consensus(base_input_path='/datasets_ieee/PDB',
    #                  target='Acetylcholinesterase',
    #                  pdb_code=('4M0E', '1YL', 604, 'A'),
    #                  base_output_path='/resultados_ieee/docking',
    #                  base_selected_mols='/datasets_ieee/ChEMBL',
    #                  mol_filename='molecules',
    #                  dock6_app_path='/home/michel/progs/dock6/',
    #                  prepare_complex=False, charge_type='am1')
    
   
    perform_consensus(base_input_path='/datasets_ieee/PDB',
                      target='Butyrylcholinesterase',
                       pdb_code=('1P0I', 'BUA', 606, 'A'),
                     base_output_path='/resultados_ieee/docking',
                      base_selected_mols='/datasets_ieee/ChEMBL',
                      mol_filename='molecules',
                      dock6_app_path='/home/michel/progs/dock6/',
                      prepare_complex=False, charge_type='am1')
    
    
    #perform_consensus(base_input_path='/datasets_ieee/PDB',
    #                  target='Beta-secretase1',
    #                  pdb_code=('5HE5', '60S', 502, 'A'),
    #                  base_output_path='/resultados_ieee/docking',
    #                  base_selected_mols='/datasets_ieee/ChEMBL',
    #                  mol_filename='molecules',
    #                  dock6_app_path='/home/michel/progs/dock6/',
    #                  prepare_complex=False, charge_type='am1')
    
   
    
