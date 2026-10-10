"""Keep MOL2 residues contiguous for DOCK6 footprint evaluation."""
from pathlib import Path


def footprint_receptor(source,destination):
    lines=Path(source).read_text().splitlines();sections={};order=[];current=None
    for line in lines:
        if line.startswith('@<TRIPOS>'):
            current=line;sections[current]=[];order.append(current)
        elif current is not None:sections[current].append(line)
    atoms=[line.split() for line in sections.get('@<TRIPOS>ATOM',[]) if line.strip()]
    if not atoms:raise ValueError('O receptor do footprint não contém átomos.')
    atoms.sort(key=lambda fields:int(fields[6]))
    renumber={int(fields[0]):index for index,fields in enumerate(atoms,1)}
    for index,fields in enumerate(atoms,1):fields[0]=str(index)
    sections['@<TRIPOS>ATOM']=[' '.join(fields) for fields in atoms]
    bonds=[]
    for line in sections.get('@<TRIPOS>BOND',[]):
        fields=line.split()
        if fields:fields[1]=str(renumber[int(fields[1])]);fields[2]=str(renumber[int(fields[2])]);bonds.append(' '.join(fields))
    sections['@<TRIPOS>BOND']=bonds
    substructures=[]
    for line in sections.get('@<TRIPOS>SUBSTRUCTURE',[]):
        fields=line.split()
        if fields:
            if int(fields[2]) in renumber:fields[2]=str(renumber[int(fields[2])])
            substructures.append(' '.join(fields))
    if '@<TRIPOS>SUBSTRUCTURE' in sections:sections['@<TRIPOS>SUBSTRUCTURE']=substructures
    text='\n'.join(line for header in order for line in [header,*sections[header]])+'\n'
    Path(destination).write_text(text)
    return Path(destination)


def footprint_rows(path):
    rows=[]
    import math
    for line in Path(path).read_text().splitlines():
        fields=line.split()
        if len(fields)==8 and fields[0]!='resname':
            try:
                numbers=[float(fields[i]) for i in (2,3,5,6)]
                if not all(math.isfinite(v) for v in numbers):raise ValueError()
                int(fields[1])
            except ValueError:raise ValueError('O footprint contém energias inválidas.') from None
            rows.append(fields)
    if not rows:raise ValueError('O DOCK6 não produziu energias por resíduo para o footprint.')
    return rows
