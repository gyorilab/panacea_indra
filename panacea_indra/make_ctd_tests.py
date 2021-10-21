import pickle
from collections import Counter
from emmaa.model_tests import StatementCheckingTest
from indra.ontology.bio import bio_ontology
from indra.statements.validate import print_validation_report
from indra.ontology.standardize import standardize_agent_name

CTD_CHEMICAL_DISEASE = '/Users/ben/data/ctd/ctd_chemical_disease.pkl'

pain = ('MESH', 'D010146')
pain_and_children = [pain] + bio_ontology.get_children(*pain)


def filter_objects(stmts, object_groundings):
    print('Filtering %d statements' % len(stmts))
    filtered_stmts = []
    for stmt in stmts:
        if set(stmt.obj.db_refs.items()) & set(object_groundings):
            filtered_stmts.append(stmt)
    print('Filtered to %d statements' % len(filtered_stmts))
    return filtered_stmts


def get_mappings():
    import gilda
    mappings = set()
    for stmt in stmts:
        subj = stmt.subj
        if 'CHEBI' in subj.db_refs:
            continue
        matches = gilda.ground(subj.name)
        chebi_match = None
        mesh_match = None
        for match in matches:
            if match.term.db == 'CHEBI':
                chebi_match = (
                match.term.db, match.term.id, match.term.entry_name)
            elif match.term.db == 'MESH':
                mesh_match = (
                match.term.db, match.term.id, match.term.entry_name)
            if chebi_match and mesh_match:
                mappings.add((chebi_match, mesh_match))
    chebi_cnts = Counter([m[0] for m in mappings])
    mesh_cnts = Counter([m[1] for m in mappings])
    with open('chebi_mesh_pred.tsv', 'w') as fh:
        for chebi, mesh in mappings:
            if chebi_cnts[chebi] == 1 and mesh_cnts[mesh] == 1:
                fh.write(
                    f'chebi\t{chebi[1]}\t{chebi[2]}\tskos:exactMatch\tmesh'
                    f'\t{mesh[1]}\t{mesh[2]}\tlexical\t0.9\t'
                    f'https://github.com/indralab/panacea_indra/blob/master/'
                    f'panacea_indra/make_ctd_tests.py\n')


if __name__ == '__main__':
    with open(CTD_CHEMICAL_DISEASE, 'rb') as fh:
        stmts = pickle.load(fh)

    pain_stmts = filter_objects(stmts, pain_and_children)
    print_validation_report(pain_stmts)
    for stmt in pain_stmts:
        for agent in stmt.real_agent_list():
            standardize_agent_name(agent, standardize_refs=True)
    with open('chemical_pain_ctd_stmts.pkl', 'wb') as fh:
        pickle.dump(pain_stmts, fh)

    pain_stmt_tests = [StatementCheckingTest(stmt) for stmt
                       in pain_stmts]
    with open('chemical_pain_ctd_tests.pkl', 'wb') as fh:
        pickle.dump(pain_stmt_tests, fh)
