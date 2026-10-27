"""Keep receptor atoms and all Vina atom-type maps in memory while screening.

Custom CLI options remain on the CLI path rather than silently losing settings.
Sessions are private to one worker and are never shared between projects.
"""
import hashlib
from pathlib import Path
import time

from .diagnostics import event, get_logger
from .ligand_cache import validate_prepared_ligand
from .progress import report_progress


class VinaSessions:
    def __init__(self):
        self.session = None
        self.key = None
        self.logger = get_logger('scientific_tools')
        self.warned = set()
        self.fingerprints = {}

    def _fallback(self, reason):
        if reason not in self.warned:
            event(self.logger, 'vina.session.fallback', reason)
            self.warned.add(reason)
        return False

    def run(self, configuration, directory):
        configuration, directory = Path(configuration), Path(directory)
        options = {}
        for line in configuration.read_text().splitlines():
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            name, separator, value = line.partition('=')
            if not separator or name.strip() in options:
                return self._fallback('Configuração personalizada: usando Vina CLI.')
            options[name.strip()] = value.strip()
        supported = {'receptor', 'ligand', 'out', 'center_x', 'center_y', 'center_z',
                     'size_x', 'size_y', 'size_z', 'exhaustiveness', 'num_modes',
                     'seed', 'verbosity', 'scoring', 'cpu', 'spacing', 'max_evals',
                     'min_rmsd', 'energy_range'}
        required = {'receptor', 'ligand', 'out', 'center_x', 'center_y', 'center_z',
                    'size_x', 'size_y', 'size_z'}
        if options.keys() - supported or not required <= options.keys() or options.get('scoring', 'vina') not in ('vina', 'vinardo'):
            return self._fallback('Opções especiais de docking: usando Vina CLI.')
        try:
            from vina import Vina
        except ImportError:
            return self._fallback('API Python do Vina indisponível: usando Vina CLI, com mapas calculados por composto.')
        def path(name):
            return (directory / options[name]).resolve()
        receptor = path('receptor')
        ligand = path('ligand')
        compound_identifier = ligand.name.removesuffix('.lig.pdbqt').removesuffix('.pdbqt')
        # Invalid floating-point geometry can make the native Vina library call
        # abort(3), bypassing Python exception handling.  Validate at the last
        # safe boundary as well as when the preparation/cache is produced.
        validate_prepared_ligand(ligand, 'pdbqt')
        center = [float(options['center_' + axis]) for axis in 'xyz']
        size = [float(options['size_' + axis]) for axis in 'xyz']
        spacing = float(options.get('spacing', .375))
        initialization = (options.get('scoring', 'vina'), int(options.get('cpu', 0)),
                          int(options.get('seed', 0)), int(options.get('verbosity', 1)),
                          float(options.get('min_rmsd', 1)))
        stat = receptor.stat()
        signature = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        previous = self.fingerprints.get(receptor)
        if previous is None or previous[0] != signature:
            previous = (signature, hashlib.sha256(receptor.read_bytes()).hexdigest())
            self.fingerprints[receptor] = previous
        key = (previous[1], tuple(center), tuple(size), spacing, initialization)
        if key != self.key:
            # Compute before loading any ligand: include every atom type, even
            # when the first molecule lacks types required by later compounds.
            event(self.logger, 'vina.maps.started', 'Calculando mapas uma vez para este receptor e sítio.', receptor_file=str(receptor))
            report_progress('Vina: calculando mapas do receptor e do sítio de docking.')
            self.session = None
            self.key = None
            session = Vina(sf_name=initialization[0], cpu=initialization[1], seed=initialization[2],
                           verbosity=initialization[3])
            session.set_receptor(str(receptor))
            session.compute_vina_maps(center=center, box_size=size, spacing=spacing)
            self.session, self.key = session, key
            event(self.logger, 'vina.maps.completed', 'Mapas prontos para reutilização pelos compostos.')
        else:
            event(self.logger, 'vina.maps.reused', 'Reutilizando receptor e mapas do sítio.')
        started = time.monotonic()
        fields = {'configuration': configuration.name, 'compound_identifier': compound_identifier}
        event(self.logger, 'vina.ligand.started', 'Docking do composto iniciado.', **fields)
        report_progress('Vina: docking de '+ligand.name+' com receptor e mapas reutilizados.')
        try:
            self.session.set_ligand_from_file(str(ligand))
            self.session.dock(exhaustiveness=int(options.get('exhaustiveness', 8)),
                              n_poses=int(options.get('num_modes', 9)),
                              min_rmsd=initialization[4], max_evals=int(options.get('max_evals', 0)))
            self.session.write_poses(str(path('out')), n_poses=int(options.get('num_modes', 9)),
                                     energy_range=float(options.get('energy_range', 3)), overwrite=True)
        except Exception:
            self.logger.exception('Falha no docking do composto.', extra={'event': 'vina.ligand.failed', **fields})
            raise
        event(self.logger, 'vina.ligand.completed', 'Docking do composto concluído.',
              duration_ms=round((time.monotonic() - started) * 1000), **fields)
        return True
