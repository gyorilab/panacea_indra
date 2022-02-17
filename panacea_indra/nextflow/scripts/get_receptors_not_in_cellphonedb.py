from indra.databases import hgnc_client
from make_ligand_receptor_database import get_genes_for_go_ids
from make_ligand_receptor_database import get_receptors
from make_ligand_receptor_database import ION_CHANNELS

if __name__ == '__main__':
    rec1 = get_receptors()
    with open('../../../go_unique_receptors.csv', 'r') as fh:
        rec2 = fh.read().splitlines()
    receptors = set(rec1) & set(rec2)
    receptors = [r for r in receptors if hgnc_client.get_hgnc_id(r)]
    olfactory = {r for r in receptors if r.startswith('OR')}
    receptors = set(receptors) - olfactory
    taste = {r for r in receptors if r.startswith('TAS')}
    receptors = set(receptors) - taste
    ion_channels = set()
    with open(ION_CHANNELS, 'r') as fh:
        for line in fh:
            ion_channels.add(line.strip())
    receptors = set(receptors) - ion_channels
    with open('receptors_not_in_cellphonedb.csv', 'w') as fh:
        fh.write('\n'.join(sorted(receptors)))
