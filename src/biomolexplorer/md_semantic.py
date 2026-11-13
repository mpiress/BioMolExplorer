"""Current-run compatibility for legacy Amber Score target descriptors."""
from pathlib import Path

def semantic_context(pipeline,results,consumer):
    """Repair missing legacy target ports in a detached current-run context.

    Mixed manifests can contain scored poses with ports and prepared target
    descriptors without them. Preserve IDs, paths and provenance; never fill
    a declared port, invent an artifact, or use another run's publications.
    """
    from .bindings import sources
    owners={ref['stage'] for ref in sources(consumer.get('bindings',{}).get('base_input_path',{}))
            if ref.get('stage') and ref.get('mode','current')=='current'}
    if not owners:return list(pipeline or [])
    from .md_inputs import receptor_input_files
    from .artifact_contracts import infer_output_port,port_types
    context=[]
    for producer in pipeline or []:
        operation=producer.get('operation')
        port=infer_output_port(operation,'receptor')
        if producer.get('id') not in owners or 'receptor_pdb' not in port_types(operation,port):
            context.append(producer);continue
        available=set(map(str,results.get(producer['id'],[])))
        records=producer.get('_execution_artifact_records') or producer.get('artifact_records') or []
        updated=[];changed=False
        for record in records:
            path=record.get('path','')
            if not record.get('output_port') and path in available and Path(path).suffix.lower()=='.pdb' and receptor_input_files([path]):
                record=dict(record,output='receptor',output_port=port,artifact_type='receptor_pdb')
                changed=True
            updated.append(record)
        context.append(dict(producer,_execution_artifact_records=updated) if changed else producer)
    return context

