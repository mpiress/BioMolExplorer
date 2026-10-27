"""Dependency-scoped checkpoints for shared DOCK6 grids and compound jobs."""
from pathlib import Path

from .calculation_cache import CalculationCache, identity
from .stage_cache import file_digest
from .paths import resolve_path


class Dock6Checkpoints:
    def __init__(self, engine, installation, receptor, settings, output):
        self.engine, self.output = engine, Path(output)
        self.receptor_identifier=receptor
        self.cache = CalculationCache('dock6', engine.export_path)
        root=Path(installation)
        self.tools={str(p.relative_to(root)):file_digest(p) for folder in ('bin','parameters')
                    for p in sorted((root/folder).iterdir()) if p.is_file()
                    and (folder=='parameters' or p.name in ('dock6','sphgen','sphere_selector','showbox','grid'))}
        self.receptor={p.name:file_digest(p) for p in sorted(Path(engine.receptorpath).glob(receptor+'.*'))
                       if p.is_file()}
        self.settings=settings
        self.templates={name:file_digest(resolve_path('src/scripts/dock6/'+name))
                        for name in ('INSPH.template','showbox.template','grid.template','min.template',
                                     'docking.template','footprint.template')}
        from .paths import SOURCE_ROOT
        self.adapter=identity({str(p.relative_to(SOURCE_ROOT)):file_digest(p) for p in
            (Path(__file__),SOURCE_ROOT/'caad/docking.py',
             Path(__file__).with_name('molecular_surface.py'),Path(__file__).with_name('footprints.py'))})
        self.hashes={}

    def fingerprint(self,path):
        stat=path.stat();stamp=(stat.st_ino,stat.st_mtime_ns,stat.st_size)
        previous=self.hashes.get(path)
        if previous is None or previous[0]!=stamp:
            previous=(stamp,file_digest(path));self.hashes[path]=previous
        return previous[1]

    def site_key(self):
        # Older direct callers select a site per ligand. Their grids depend on
        # those coordinates; the platform uses a fixed reference center.
        ligands={} if self.settings['center'] is not None else {
            p.name:self.fingerprint(p) for p in Path(self.engine.ligandpath).glob('*.mol2')}
        return identity(dict(receptor=self.receptor,settings=self.settings,tools=self.tools,adapter=self.adapter,
                             ligands=ligands,
                             templates={name:self.templates[name] for name in ('INSPH.template','showbox.template','grid.template')}))

    def site_files(self):
        return {p.relative_to(self.output).as_posix():p
                for directory in ('surface','showbox','gridbox')
                for p in (self.output/directory).rglob('*')
                if p.is_file() and (p.suffix in ('.sph','.dms','.nrg','.bmp') or p.name.endswith('.box.pdb'))}

    def restore_site(self):
        key=self.site_key();entry=self.cache.metadata(key)
        if not entry:return False
        if any(not (self.output/name).resolve().is_relative_to(self.output.resolve()) for name in entry.get('files',[])):
            return False
        return self.cache.restore(key,{name:self.output/name for name in entry['files']}) is not None

    def remember_site(self):
        files=self.site_files()
        if files:self.cache.remember(self.site_key(),files,dict(files=list(files),kind='receptor-site',
                                                              receptor_identifier=self.receptor_identifier))

    def compound_key(self, method, code):
        template={'perform_parallel_minimization':'min.template',
                  'perform_parallel_docking':'docking.template',
                  'perform_parallel_footprint':'footprint.template'}[method]
        # Only this compound's dependencies enter the key. Adding a different
        # ligand does not change a completed minimization/docking/footprint.
        dependencies={}
        for directory in (Path(self.engine.ligandpath),self.output/'energy_min',self.output/'flex',self.output/'rigid'):
            if method=='perform_parallel_minimization' and directory!=Path(self.engine.ligandpath):continue
            if method=='perform_parallel_docking' and directory in (self.output/'flex',self.output/'rigid'):continue
            for p in sorted(directory.glob(code+'*.mol2')):
                from .dock6_validation import compound_code
                if compound_code(p)==code:dependencies[str(directory.name)+'/'+p.name]=self.fingerprint(p)
        grid={p.relative_to(self.output).as_posix():self.fingerprint(p) for p in self.site_files().values()}
        return identity(dict(site=self.site_key(),grid=grid,method=method,compound=code,inputs=dependencies,
                             template=self.templates[template],mode=self.engine._Dock6__conformer_search_type
                             if method=='perform_parallel_docking' else None))
