"""Topology-defined, streaming analysis of an implicit-solvent Amber trajectory."""
from pathlib import Path
import re

import numpy as np
from scipy.spatial import cKDTree

from .storage import write_json

PROTEIN=set('ALA ARG ASN ASP ASH CYS CYM CYX GLN GLU GLH GLY HIS HID HIE HIP ILE LEU LYS LYN MET PHE PRO SER THR TRP TYR VAL'.split())


def topology(path):
    """Read Amber's fixed-width topology; do not infer atoms from filenames."""
    sections={};name=None;width=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('%FLAG '):name=line.split()[1];sections[name]=[];width=None
        elif line.startswith('%FORMAT'):
            match=re.search(r'\(\d+([aAiIeEfF])(\d+)',line)
            if not match:raise ValueError('Unsupported Amber topology format: '+line)
            kind,width=match[1].lower(),int(match[2])
        elif name and width and not line.startswith('%'):
            for start in range(0,len(line),width):
                cell=line[start:start+width].strip()
                if cell:sections[name].append(cell if kind=='a' else int(cell) if kind=='i' else float(cell.replace('D','E')))
    pointers=sections.get('POINTERS',[])
    if len(pointers)<12:raise ValueError('Incomplete Amber topology POINTERS.')
    count,nres=pointers[0],pointers[11]
    for flag in ('ATOM_NAME','MASS','CHARGE'):
        if len(sections.get(flag,[]))!=count:raise ValueError('Incomplete Amber topology: '+flag)
    if len(sections.get('RESIDUE_POINTER',[]))!=nres or len(sections.get('RESIDUE_LABEL',[]))!=nres:
        raise ValueError('Incomplete Amber residue map.')
    masses=np.asarray(sections['MASS'],float);charges=np.asarray(sections['CHARGE'],float)/18.2223
    if not np.isfinite(masses).all() or not np.isfinite(charges).all() or (masses<=0).any():
        raise ValueError('Invalid masses or charges in Amber topology.')
    residue=np.empty(count,dtype=int)
    ends=sections['RESIDUE_POINTER'][1:]+[count+1]
    for index,(start,end) in enumerate(zip(sections['RESIDUE_POINTER'],ends)):
        if not 1<=start<end<=count+1:raise ValueError('Invalid Amber residue boundaries.')
        residue[start-1:end-1]=index
    return dict(natom=count,nres=nres,names=sections['ATOM_NAME'],masses=masses,charges=charges,
                residue=residue,labels=sections['RESIDUE_LABEL'],
                bonds=[(int(values[i])//3,int(values[i+1])//3)
                       for values in (sections.get('BONDS_INC_HYDROGEN',[]),sections.get('BONDS_WITHOUT_HYDROGEN',[]))
                       for i in range(0,len(values),3)])


def selection(directory,require_protein=True):
    directory=Path(directory)
    total=topology(directory/'rec.lig.1.prmtop');lig=topology(directory/'lig.1.prmtop');rec=topology(directory/'rec.prmtop')
    for prefix,info in [('rec.lig.1',total),('lig.1',lig),('rec',rec)]:
        lines=(directory/(prefix+'.inpcrd')).read_text().splitlines()
        try:
            if int(lines[1].split()[0])!=info['natom']:raise ValueError()
            values=[float(line[i:i+12]) for line in lines[2:2+(info['natom']*3+5)//6] for i in range(0,len(line),12) if line[i:i+12].strip()]
            if len(values)!=info['natom']*3 or not np.isfinite(values).all():raise ValueError()
        except (ValueError,IndexError):raise ValueError('Invalid Amber coordinates for '+prefix) from None
    n=lig['natom']
    if n+rec['natom']!=total['natom'] or total['names'][:n]!=lig['names'] or total['names'][n:]!=rec['names']:
        raise ValueError('Amber complex topology does not match ligand-first/receptor atom order.')
    heavy=total['masses']>2.5
    ligand=np.flatnonzero(heavy & (np.arange(total['natom'])<n))
    protein=np.flatnonzero(heavy & (np.arange(total['natom'])>=n) &
        np.array([total['labels'][r] in PROTEIN for r in total['residue']]))
    backbone=np.array([i for i in protein if total['names'][i] in ('N','CA','C','O')],dtype=int)
    ca=np.array([i for i in protein if total['names'][i]=='CA'],dtype=int)
    if len(ligand)==0 or require_protein and (len(backbone)<3 or len(ca)==0):
        raise ValueError('Trajectory analysis requires ligand heavy atoms and a protein backbone with C-alpha atoms.')
    return total,dict(ligand=ligand,protein=protein,backbone=backbone,ca=ca,
                      ligand_mask=f':1-{lig["nres"]}',protein_mask=':'+','.join(str(r+1) for r in sorted(set(total['residue'][protein]))))


def fit_frame(frame,reference,indices):
    moving=frame[indices];target=reference[indices]
    center=moving.mean(axis=0);target_center=target.mean(axis=0)
    u,_,vt=np.linalg.svd((moving-center).T@(target-target_center))
    rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    return (frame-center)@rotation+target_center


def radius(coords,masses):
    center=np.average(coords,axis=0,weights=masses)
    return float(np.sqrt(np.average(np.sum((coords-center)**2,axis=1),weights=masses)))


def frames(path):
    from scipy.io import netcdf_file
    # Amber's ioutfm=1 trajectory is classic NetCDF, which SciPy can stream
    # without loading the complete trajectory into RAM.
    with netcdf_file(path,'r',mmap=True) as stream:
        variable=stream.variables['coordinates']
        for index in range(variable.shape[0]):yield np.array(variable[index],dtype=float,copy=True)
        del variable


def analyze_frames(sequence,info,masks,interval_ps):
    first=None;series=[];mean=None;m2=None;occupancy={};site=None;count=0
    for count,frame in enumerate(sequence,1):
        frame=np.asarray(frame,float)
        if frame.shape!=(info['natom'],3) or not np.isfinite(frame).all():raise ValueError('Invalid or non-finite trajectory coordinates.')
        if first is None:
            first=frame.copy()
            near=cKDTree(first[masks['ligand']]).query(first[masks['protein']])[0]<=5.0
            site_residues=set(info['residue'][masks['protein'][near]])
            site=np.array([i for i in masks['protein'] if info['residue'][i] in site_residues],dtype=int)
            mean=np.zeros_like(frame);m2=np.zeros_like(frame)
        fitted=fit_frame(frame,first,masks['backbone'])
        delta=fitted-mean;mean+=delta/count;m2+=delta*(fitted-mean)
        rms=lambda mask:float(np.sqrt(np.mean(np.sum((fitted[mask]-first[mask])**2,axis=1))))
        row=dict(frame=count,time_ps=count*interval_ps,protein_backbone_rmsd=rms(masks['backbone']),
                 ligand_rmsd=rms(masks['ligand']),protein_rmsd=rms(masks['protein']),
                 complex_rg=radius(frame,info['masses']),
                 protein_rg=radius(frame[masks['protein']],info['masses'][masks['protein']]))
        if len(site):
            ligand_center=np.average(frame[masks['ligand']],axis=0,weights=info['masses'][masks['ligand']])
            site_center=np.average(frame[site],axis=0,weights=info['masses'][site])
            row['ligand_site_distance']=float(np.linalg.norm(ligand_center-site_center))
        distances=cKDTree(frame[masks['ligand']]).query(frame[masks['protein']])[0]
        for residue in set(info['residue'][masks['protein'][distances<=4.0]]):
            occupancy[int(residue)]=occupancy.get(int(residue),0)+1
        series.append(row)
    if count<2:raise ValueError('At least two saved trajectory frames are required for fluctuation analysis; increase steps or reduce stride.')
    fluctuations=np.sqrt(np.sum(m2/count,axis=1))
    rmsf=[dict(residue=int(info['residue'][i])+1,name=info['labels'][info['residue'][i]],
               rmsf=float(fluctuations[i]),binding_site=bool(len(site) and info['residue'][i] in set(info['residue'][site]))) for i in masks['ca']]
    contacts=[dict(residue=r+1,name=info['labels'][r],occupancy=count_value/count) for r,count_value in sorted(occupancy.items())]
    summary=dict(frames=count,sampled_duration_ps=(count-1)*interval_ps,first_frame_ps=interval_ps,
                 last_frame_ps=count*interval_ps,
                 most_flexible_residues=sorted(rmsf,key=lambda r:r['rmsf'],reverse=True)[:10],
                 most_persistent_contacts=sorted(contacts,key=lambda r:r['occupancy'],reverse=True)[:10])
    for name in series[0]:
        if name in ('frame','time_ps'):continue
        summary['mean_'+name]=float(np.mean([r[name] for r in series]));summary['final_'+name]=series[-1][name]
    return dict(version=1,summary=summary,series=series,rmsf=rmsf,contacts=contacts,
        definitions=dict(reference='First saved frame; unweighted protein-backbone rigid fit',
            ligand_rmsd='Ligand heavy atoms after protein fit, without ligand refitting',
            rmsf='Protein C-alpha RMSF about the trajectory mean after protein-backbone fitting',
            complex_rg='Mass-weighted radius of gyration, all complex atoms',
            protein_rg='Mass-weighted protein-heavy-atom radius of gyration',
            contacts='Fraction of saved frames with any protein-residue/ligand heavy-atom distance <=4 Angstrom',
            ligand_site_distance='Mass-weighted ligand/site COM distance; static site consists of protein residues within 5 Angstrom in the first saved frame'),
        masks={name:value.tolist() if isinstance(value,np.ndarray) else value for name,value in masks.items()})



def coordinates(path,count):
    lines=Path(path).read_text().splitlines()
    try:
        if int(lines[1].split()[0])!=count:raise ValueError()
        values=[float(line[i:i+12]) for line in lines[2:2+(3*count+5)//6] for i in range(0,len(line),12) if line[i:i+12].strip()]
        result=np.asarray(values).reshape(count,3)
        if not np.isfinite(result).all():raise ValueError()
    except (ValueError,IndexError):raise ValueError('Invalid Amber coordinate count or non-finite coordinates: '+Path(path).name) from None
    return result


def analyze_final_state(directory):
    from .md_engine import plot
    directory=Path(directory);info,masks=selection(directory,require_protein=False)
    initial=coordinates(directory/'rec.lig.1.inpcrd',info['natom'])
    final_path=directory/'rec.lig.1.final_pose.amber.restart'
    if not final_path.is_file():
        data=dict(version=2,summary=dict(initial_complex_rg=radius(initial,info['masses']),hydrogen_atoms=int(np.count_nonzero(info['masses']<2.5))),
            definitions=dict(initial_complex_rg='Mass-weighted radius of gyration of the prepared initial complex only.',hydrogen_atoms='Hydrogen/isotope atom count inferred from topology masses below 2.5 u; not hydrogen bonds.'),
            unavailable='Native DOCK6 did not write a final complex; structural comparison is unavailable. No trajectory or RMSF is produced.')
    else:
        final=coordinates(final_path,info['natom'])
        fitted=fit_frame(final,initial,masks['backbone']) if len(masks['backbone'])>=3 else None
        rms=lambda mask:float(np.sqrt(np.mean(np.sum((fitted[mask]-initial[mask])**2,axis=1))))
        summary=dict(initial_complex_rg=radius(initial,info['masses']),final_complex_rg=radius(final,info['masses']),hydrogen_atoms=int(np.count_nonzero(info['masses']<2.5)))
        if fitted is not None:summary.update(final_protein_backbone_rmsd=rms(masks['backbone']),final_ligand_rmsd=rms(masks['ligand']))
        data=dict(version=2,summary=summary,definitions=dict(final_ligand_rmsd='Heavy atoms versus prepared initial complex after protein-backbone fit; no ligand refitting or symmetry correction.',final_protein_backbone_rmsd='Protein N/CA/C/O atoms versus prepared initial complex after rigid backbone fit.',complex_rg='Mass-weighted radius of gyration of all complex atoms, initial/final structures only.',hydrogen_atoms='Hydrogen/isotope atom count inferred from topology masses below 2.5 u; not hydrogen bonds.'),unavailable='No native trajectory: RMSF, time series, hydrogen-bond occupancy and equilibrium claims are unavailable.')
        if fitted is not None:plot(directory/'structural_rmsd.png',[0,1],[summary['final_protein_backbone_rmsd'],summary['final_ligand_rmsd']],'Initial versus final structure RMSD','0: backbone · 1: ligand','Angstrom')
        plot(directory/'structural_rg.png',[0,1],[summary['initial_complex_rg'],summary['final_complex_rg']],'Complex radius of gyration','0: initial · 1: final','Angstrom')
    write_json(data,directory/'analysis.json');return data
