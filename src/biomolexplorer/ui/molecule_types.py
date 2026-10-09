"""ChEMBL molecule_type values, stored in the crawler's lowercase format.

Verified against https://www.ebi.ac.uk/chembl/api/data/molecule.json
using molecule_type__iexact and only=molecule_type (October 2026).
Keep the vocabulary local so opening a form does not require network access.
"""

MOLECULE_TYPES = (
    'Small molecule', 'Antibody', 'Antibody drug conjugate', 'Cell', 'Enzyme',
    'Gene', 'Oligonucleotide', 'Oligosaccharide', 'Protein', 'Vaccine component',
    'Unknown',
)
