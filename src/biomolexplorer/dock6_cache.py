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
        import ast
        tree=ast.parse((SOURCE_ROOT/'caad/docking.py').read_text())
        cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='Dock6')
        methods={node.name:ast.dump(node,include_attributes=False) for node in cls.body
                 if isinstance(node,ast.FunctionDef) and node.name in ('prepare_surface','prepare_showbox','prepare_gridbox','_grid_prefix')}
        self.site_adapter=identity(dict(surface=file_digest(Path(__file__).with_name('molecular_surface.py')),methods=methods))
        self.hashes={}

    def fingerprint(self,path):
        stat=path.stat();stamp=(stat.st_ino,stat.st_mtime_ns,stat.st_size)
        previous=self.hashes.get(path)
        if previous is None or previous[0]!=stamp:
            previous=(stamp,file_digest(path));self.hashes[path]=previous
        return previous[1]

    def site_request(self,adapter,*,legacy=False):
        # Older direct callers select a site per ligand. Their grids depend on
        # those coordinates; the platform uses a fixed reference center.
        ligands={} if self.settings['center'] is not None else {
            p.name:self.fingerprint(p) for p in Path(self.engine.ligandpath).glob('*.mol2')}
        receptor=self.receptor if legacy else {name:value for name,value in self.receptor.items()
            if name.endswith(('.dockprep.mol2','.noH.pdb'))}
        return identity(dict(receptor=receptor,settings=self.settings,tools=self.tools,adapter=adapter,
                             ligands=ligands,
                             templates={name:self.templates[name] for name in ('INSPH.template','showbox.template','grid.template')}))

    def site_key(self):
        return self.site_request(self.site_adapter)

    def site_files(self):
        return {p.relative_to(self.output).as_posix():p
                for directory in ('surface','showbox','gridbox')
                for p in (self.output/directory).rglob('*')
                if p.is_file() and (p.suffix in ('.sph','.dms','.nrg','.bmp') or p.name.endswith('.box.pdb'))}

    def restore_site(self):
        import json
        keys=[self.site_key(),self.site_request(self.site_adapter,legacy=True)]
        if self.settings['center'] is not None:
            from .templates import RESOURCE_ROOT
            compatibility=json.loads((RESOURCE_ROOT/'dock6_cache_compatibility.json').read_text())
            keys.extend(self.site_request(adapter,legacy=True) for adapter in compatibility['legacy_site_adapters'])
        for key in keys:
            entry=self.cache.metadata(key)
            if not entry or not entry.get('files'):continue
            if any(not (self.output/name).resolve().is_relative_to(self.output.resolve()) for name in entry['files']):continue
            if self.cache.restore(key,{name:self.output/name for name in entry['files']}) is not None:
                if key!=keys[0]:self.remember_site()
                return True
        return False

    def remember_site(self):
        files=self.site_files()
        if files:self.cache.remember(self.site_key(),files,dict(files=list(files),kind='receptor-site',
                                                              receptor_identifier=self.receptor_identifier))

    def compound_key(self, method, code):
        template={'perform_parallel_minimization':'min.template',
                  'perform_parallel_docking':'docking.template',
                  'perform_parallel_footprint':'footprint.template'}[method]
        mode=self.engine.protocol_for(code)
        if method=='perform_parallel_docking' and mode not in ('denovo','redock'):template='min.template'
        # Only this compound's dependencies enter the key. Adding a different
        # ligand does not change a completed minimization/docking/footprint.
        dependencies={}
        for directory in (Path(self.engine.ligandpath),self.output/'energy_min',self.output/'flex',self.output/'rigid'):
            if method=='perform_parallel_minimization' and directory!=Path(self.engine.ligandpath):continue
            if method=='perform_parallel_docking' and directory in (self.output/'flex',self.output/'rigid'):continue
            if method=='perform_parallel_docking' and mode not in ('denovo','redock') and directory!=Path(self.engine.ligandpath):continue
            for p in sorted(directory.glob(code+'*.mol2')):
                from .dock6_validation import compound_code
                if compound_code(p)==code:dependencies[str(directory.name)+'/'+p.name]=self.fingerprint(p)
        grid={p.relative_to(self.output).as_posix():self.fingerprint(p) for p in self.site_files().values()}
        from .dock6_protocol import policy
        protocol=dict(mode=mode,parameters=policy(mode)) if method=='perform_parallel_docking' else None
        return identity(dict(site=self.site_key(),grid=grid,method=method,compound=code,inputs=dependencies,
                             template=self.templates[template],protocol=protocol,
                             mode=self.engine._Dock6__conformer_search_type
                             if method=='perform_parallel_docking' and mode in ('denovo','redock') else None))
