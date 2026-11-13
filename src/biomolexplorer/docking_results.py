"""Authorized score tables with audited row removal and actual-pose previews."""
from pathlib import Path
from .compound_tables import CompoundTables
from .input_validation import columns
from .workspace import AccessDenied


class DockingResults(CompoundTables):
    available_statuses=('succeeded','running','failed','cancelled','interrupted')
    def _allowed(self,pid,stage):
        key=(pid,self._stage_snapshot[0] if self._stage_snapshot else None,stage.get('id'),tuple(stage.get('artifacts',[])))
        cached=getattr(self,'_allowed_snapshot',None)
        if cached and cached[0]==key:return cached[1]
        allowed={str(self.store.scoped_path(pid,name)) for name in stage.get('artifacts',[])}
        self._allowed_snapshot=(key,allowed)
        return allowed

    def page(self,token,pid,rid,sid,filename,*args,**kwargs):
        result=super().page(token,pid,rid,sid,filename,*args,**kwargs)
        path=self._path(token,pid,rid,sid,filename)
        _,stage=self._stage(token,pid,rid,sid)
        if stage.get('operation')=='docking_vina':return result
        allowed=self._allowed(pid,stage)
        for row in result['rows']:
            try:
                self._footprint_path(pid,path,row['properties'],stage,allowed=allowed)
                row['footprint_available']=True
            except (ValueError,AccessDenied):row['footprint_available']=False
        return result

    def tables(self,token,pid,rid,sid):
        _,stage=self._stage(token,pid,rid,sid)
        result=[]
        artifacts=stage.get('artifacts',[]) if stage['status']=='succeeded' else list(dict.fromkeys(
            path for batch in stage.get('batches',[]) if batch.get('status')=='succeeded'
            for path in batch.get('artifacts',[])))
        for name in artifacts:
            if not str(name).lower().endswith('.csv'):continue
            path=self.store.scoped_path(pid,name)
            if path.suffix=='.csv' and path.is_file() and {'molecule_chembl_id','canonical_smiles'}<=columns(path):
                fields=columns(path)
                if 'score' in fields or {'vina','dock6'}<=fields:
                    result.append(dict(path=str(path),name=path.name,integrated=False))
        # Keep a summary for every batch, without duplicated receptor tables.
        from .docking_data import result_tables
        authoritative={str(p) for p in result_tables(t['path'] for t in result)}
        result=[t for t in result if t['path'] in authoritative]
        from .artifacts import annotate_files
        return annotate_files(stage,result,rid)

    def _target_route(self,token,pid,rid,sid,table,index,version,column,target):
        from .consensus_targets import target_route
        return target_route(self,token,pid,rid,sid,table,index,version,target,column)

    def target_choices(self,token,pid,rid,sid,table,index,version,column='vina_pose',action='pose',validate=True):
        from .consensus_targets import target_details,target_label
        if column not in ('vina_pose','dock6_pose','conformer_file') or action not in ('pose','footprint'):
            raise AccessDenied('Conformação não autorizada.')
        if type(index) is not int or index<0:raise AccessDenied('Conformação não autorizada.')
        path=self._path(token,pid,rid,sid,table)
        row=self._record(path,index,version,root=self.store.project_dir(pid))
        run,stage=self._stage(token,pid,rid,sid)
        result=[]
        for target in target_details(row):
            available=True
            try:
                if not validate:available=True
                elif action=='footprint':self.footprint(token,pid,rid,sid,table,index,version,target=target)
                else:self.scene(token,pid,rid,sid,table,index,version,column=column,target=target)
            except (OSError,ValueError,AccessDenied):available=False
            result.append(dict(id=target,label=target_label(self,pid,run,row,target),available=available))
        self.store.project(token,pid)
        return result

    def pose(self,token,pid,rid,sid,table,index,version,column='conformer_file',target=None):
        if column not in ('conformer_file','vina_pose','dock6_pose'):raise AccessDenied('Conformação não autorizada.')
        if target is not None:
            route=self._target_route(token,pid,rid,sid,table,index,version,column,target)
            if route:
                child,source_rid,source_sid,selection=route
                return child.pose(token,pid,source_rid,source_sid,**selection)
        path=self._path(token,pid,rid,sid,table)
        if type(index) is not int or index<0 or column not in ('conformer_file','vina_pose','dock6_pose'):raise AccessDenied('Conformação não autorizada.')
        try:row=self._record(path,index,version,root=self.store.project_dir(pid))
        except IndexError:raise AccessDenied('Conformação não autorizada.') from None
        value=row.get(column)
        if not value:raise ValueError('A conformação selecionada não está disponível.')
        pose=self.store.scoped_path(pid,Path(value) if Path(value).is_absolute() else path.parent/value)
        _,stage=self._stage(token,pid,rid,sid)
        if str(pose) not in stage.get('artifacts',[]) or not pose.is_file():raise AccessDenied('Conformação não autorizada para esta etapa.')
        return str(pose)

    def _row(self,token,pid,rid,sid,table,index,version,column):
        pose=Path(self.pose(token,pid,rid,sid,table,index,version,column))
        path=self._path(token,pid,rid,sid,table)
        try:row=self._record(path,index,version,root=self.store.project_dir(pid))
        except IndexError:raise AccessDenied('Conformação não autorizada.') from None
        _,stage=self._stage(token,pid,rid,sid)
        return pose,path,row,stage

    def _context(self,token,pid,rid,sid,table,index,version,column):
        pose,path,row,stage=self._row(token,pid,rid,sid,table,index,version,column)
        artifacts=self._allowed(pid,stage)
        batches=[b for b in stage.get('batches',[]) if str(path) in b.get('artifacts',[])]
        if len(batches)==1:contexts=[batches[0]]
        elif not batches:contexts=[stage]
        else:raise AccessDenied('Conformação não autorizada.')
        inputs=set()
        for context in contexts:
            for name in context.get('input_files',{}).get('base_input_path',[]):
                inputs.add(str(self.store.scoped_path(pid,name)))
        allowed=artifacts|inputs
        def linked(field):
            value=row.get(field)
            if not value:return None
            candidate=self.store.scoped_path(pid,Path(value) if Path(value).is_absolute() else path.parent/value)
            if str(candidate) not in allowed or not candidate.is_file():raise AccessDenied('Conformação não autorizada para esta etapa.')
            return candidate
        engine='vina' if column=='vina_pose' else 'dock6' if column=='dock6_pose' else row.get('engine','')
        receptor=linked(engine+'_receptor_file') or linked('receptor_file')
        if receptor is None:
            receptor_id=row.get('receptor_id','')
            candidates=[self.store.scoped_path(pid,name) for name in allowed if Path(name).name in
                (receptor_id+'.noH.pdb',receptor_id+'.dockprep.pdbqt',receptor_id+'.dockprep.mol2')]
            for suffix in ('.noH.pdb','.dockprep.pdbqt','.dockprep.mol2'):
                selected=[p for p in candidates if p.name.endswith(suffix) and p.is_file()]
                if len(selected)==1:receptor=selected[0];break
                if len(selected)>1:raise ValueError('O receptor associado a esta pose é ambíguo.')
        return pose,path,row,stage,receptor,linked,engine

    def scene(self,token,pid,rid,sid,table,index,version,column='conformer_file',target=None):
        if column not in ('conformer_file','vina_pose','dock6_pose'):raise AccessDenied('Conformação não autorizada.')
        if target is not None:
            route=self._target_route(token,pid,rid,sid,table,index,version,column,target)
            if route:
                child,source_rid,source_sid,selection=route
                return child.scene(token,pid,source_rid,source_sid,**selection)
        pose,path,row,stage,receptor,linked,engine=self._context(token,pid,rid,sid,table,index,version,column)
        if receptor is None:raise ValueError('O receptor associado a esta pose não está disponível.')
        reference=linked(engine+'_reference_file') or linked('reference_file')
        layers=[dict(role='receptor',path=str(receptor),format=receptor.suffix[1:]),
            dict(role='pose',path=str(pose),format=pose.suffix[1:])]
        if reference:layers.append(dict(role='reference',path=str(reference),format=reference.suffix[1:]))
        return dict(name=row['molecule_chembl_id']+' · '+row.get('receptor_id','')+' · '+engine.upper(),
            layers=layers,reference_kind='prepared' if reference else None,ligand_smiles=row.get('canonical_smiles'))

    def footprint(self,token,pid,rid,sid,table,index,version,target=None):
        if target is not None:
            route=self._target_route(token,pid,rid,sid,table,index,version,'footprint_file',target)
            if route:
                child,source_rid,source_sid,selection=route
                selection.pop('column')
                return child.footprint(token,pid,source_rid,source_sid,**selection)
        path=self._path(token,pid,rid,sid,table)
        if type(index) is not int or index<0:raise AccessDenied('Footprint não autorizado.')
        try:row=self._record(path,index,version,root=self.store.project_dir(pid))
        except IndexError:raise AccessDenied('Footprint não autorizado.') from None
        _,stage=self._stage(token,pid,rid,sid)
        return self._footprint_path(pid,path,row,stage)

    def footprint_pdf(self,token,pid,rid,sid,table,index,version,target=None):
        from .visualizations import MAX_VIEW_BYTES
        path=self.footprint(token,pid,rid,sid,table,index,version,target=target)
        data=self.store.read_file(token,pid,path,MAX_VIEW_BYTES+1)
        if len(data)>MAX_VIEW_BYTES:raise ValueError('O PDF ultrapassa o limite de 32 MB.')
        return data

    def _footprint_path(self,pid,path,row,stage,allowed=None):
        if allowed is None:allowed=self._allowed(pid,stage)
        def authorized(candidate):
            # Native runs historically did not list PDFs in their manifests.
            native=stage.get('operation') in ('docking_dock6','consensus') and not stage.get('configuration',{}).get('provided_results')
            return str(candidate) in allowed or native and candidate.is_relative_to(path.parent.resolve())
        def checked(candidate):
            companion=candidate.parent.parent/(row['molecule_chembl_id']+'_footprint_scored.txt')
            companion=self.store.scoped_path(pid,companion)
            if authorized(companion) and companion.is_file():
                from .footprints import footprint_rows
                from .result_index import file_stamp
                stamp=(str(companion),file_stamp(companion))
                validated=getattr(self,'_validated_footprints',{})
                if stamp not in validated:
                    footprint_rows(companion)
                    if len(validated)>=128:validated.pop(next(iter(validated)))
                    validated[stamp]=True;self._validated_footprints=validated
            return str(candidate)
        value=row.get('footprint_file')
        if value:
            candidate=self.store.scoped_path(pid,Path(value) if Path(value).is_absolute() else path.parent/value)
            if candidate.suffix.lower()=='.pdf' and authorized(candidate) and candidate.is_file():return checked(candidate)
            raise AccessDenied('Footprint não autorizado para esta etapa.')
        # Earlier native DOCK6 runs store plots next to the rigid/flex pose folder.
        value=row.get('dock6_pose') or row.get('conformer_file')
        if not value:raise ValueError('O gráfico footprint associado a esta pose não está disponível.')
        pose=self.store.scoped_path(pid,Path(value) if Path(value).is_absolute() else path.parent/value)
        expected=self.store.scoped_path(pid,pose.parent.parent/'footprint/plots'/(row['molecule_chembl_id']+'.pdf'))
        if authorized(expected) and expected.is_file():return checked(expected)
        raise ValueError('O gráfico footprint associado a esta pose não está disponível.')
