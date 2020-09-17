import pickle
from indra.assemblers.cx import CxAssembler
from indra.assemblers.html import HtmlAssembler
from indra.assemblers.cx.hub_layout import add_semantic_hub_layout


def assemble_html(stmts_with_counts):
    for channel, (stmts, ev_counts, source_counts) in stmts_with_counts.items():
        ha = HtmlAssembler(stmts, ev_totals=ev_counts,
                           source_counts=source_counts,
                           title='INDRA statements for %s' % channel,
                           db_rest_url='https://db.indra.bio')
        ha.save_model('../html/%s.html' % channel)


def assemble_cx(stmts_with_counts):
    for channel, (stmts, ev_counts, source_counts) in stmts_with_counts.items():
        cxa = CxAssembler(stmts, network_name='%s INDRA interactome' % channel)
        cxa.make_model()
        add_semantic_hub_layout(cxa.cx, channel)
        model_id = cxa.upload_model()


if __name__ == '__main__':
    with open('../ion_channel_stmts_v4.pkl', 'rb') as fh:
        stmts_with_counts = pickle.load(fh)
    for channel, (stmts, ev_counts, source_counts) in stmts_with_counts.items():
        for stmt in stmts:
            ev_counts[stmt.get_hash()] = len(stmt.evidence)
    #assemble_cx(stmts_with_counts)
    assemble_html(stmts_with_counts)
