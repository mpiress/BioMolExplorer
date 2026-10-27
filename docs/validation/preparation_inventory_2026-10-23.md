# Inventário da preparação molecular (antes da migração)

## Caminhos inspecionados

- `docking_preparation._prepare_one`: moléculas externas entram por SMILES ou pose; Open Babel escreve MOL2, Chimera executa remoção/adição de H, `addcharge ... chargeModel 14sb method gas|am1` e minimização opcional; Open Babel exporta PDBQT. Um Chimera é iniciado por ligante não reutilizado. `LigandCache` preserva produtos por composto/opções/formato.
- `caad.docking.Docking.prepare_for_docking`: scripts Chimera separam complexo/cadeia, receptor/cofatores e ligante cristalográfico. `prepare_on_obabel` exporta PDBQT rígido do receptor. Nomes `.complex.pdb`, `.dockprep.mol2`, `.dockprep.pdbqt`, `.noH.pdb`, centros e `pdb_codes.csv` constituem contratos existentes. `StructureCache` reutiliza receptor por PDB/cadeia/opções e ligantes por par.
- `caad.redocking` contém o caminho CLI legado equivalente. `wrappers.redocking` utiliza os checkpoints de preparação e resultados incrementais.
- `docking_data.prepare_ligands`: consome produtos preparados ou converte poses/SMILES. Há um caminho histórico adicional de AM1 via Chimera e Gasteiger via Open Babel, separado do preparador de candidatos.
- Templates em `resources/chimera`: `prepare_complex`, `prepare_receptor`, `prepare_ligand` e `prepare_better_conform`. O receptor usa templates AMBER 14SB para resíduos padrão; `gas`/`am1` controlam resíduos não padrão, **não** uma carga Gasteiger/AM1 da proteína inteira. Cofatores selecionados são mantidos no receptor; águas/H dependem das opções por papel.
- `vina_session.VinaSessions`: API Python do Vina mantém receptor/mapas por sítio; `vina_batch` executa com isolamento contra abortos nativos. Configurações especiais continuam no CLI. Entradas/saídas PDBQT são artefatos auditáveis.
- `wrappers.docking`, `dock6_cache` e `caad.docking`: MOL2, superfície SES nativa, sphgen, seleção de esferas, showbox, grid e dock6. Preparação dos mapas é independente dos compostos. Não alterar esses parâmetros para ganhar velocidade.
- `md_engine`: `prepare_amber.pl` e helpers instalados executam antechamber, parmchk, tleap/teLeap para receptor, ligante e complexo. O bundle DOCK6 testado usa ff94/GAFF. DOCK6 AMBER Score executa min–MD–min dos três sistemas; uma trajetória sander e análises cpptraj são opcionais e distintas. Não há chamada a pmemd nem mudança implícita para GAFF2.
- `md_cache` e `calculation_cache`: checkpoints compartilhados indexados por SHA-256/SQLite; objetos imutáveis e cópias independentes/reflinks. `md_analysis` valida topologia/coordenadas e analisa trajetória. Manter o protocolo científico existente.

## pH, temporários e riscos

O pH é aplicado atualmente pelo Open Babel antes/depois de algumas exportações; Chimera AddH faz escolhas próprias. Reaplicar protonação após cargas é um risco de inconsistência e deve ser evitado no novo caminho. Não é correto alegar que o pH sozinho define todos os microestados. O legado mantém scripts `.com`, MOL2/PDBQT científicos e metadados; diretórios de candidatos são temporários. Amber exige arquivos nativos, auxiliares e logs; seus checkpoints permitem reutilização sem repetir parametrização.

Não substituir minimização AMBER do Chimera por MMFF/UFF silenciosamente. Não aplicar Gasteiger genérico a metais, centros coordenados ou proteínas. Novas implementações devem falhar explicitamente nesses casos ou manter o backend legado selecionado. Não remover águas/cofatores para fazer uma estrutura passar na validação.

## Fontes metodológicas

- [Chimera Add Charge](https://www.cgl.ucsf.edu/chimera/1.5.1/docs/ContributedSoftware/addcharge/addcharge.html): distinção entre templates de resíduos padrão e Antechamber para outros resíduos.
- [Meeko API de preparação](https://apidocmeeko.readthedocs.io/en/latest/_modules/meeko/preparation.html): `charge_model="read"` preserva cargas fornecidas, separando carga de tipagem AutoDock.
- [RDKit Gasteiger](https://rdkit.org/docs/source/rdkit.Chem.rdchem.html): cargas de átomos e hidrogênios implícitos exigem tratamento consistente.

## Estado da migração

Manter `legacy_chimera` para projetos e protocolos atuais até existir comparação representativa de receptor, docking e redocking. O novo núcleo será disponibilizado explicitamente; não há autorização científica para substituir minimização ou estados de protonação apenas por uma medição de velocidade.
