import sys
from fnvhash import fnv1a_32
from flask_wtf import FlaskForm
from flask_bootstrap import Bootstrap
from flask import Flask, render_template, request, redirect, session
from wtforms import SubmitField, validators, SelectMultipleField
from panacea_indra.get_channel_modulators import find_reg_non_regs, \
    get_agent_urls

app = Flask(__name__)
app.config['SECRET_KEY'] = 'dev_key'

Bootstrap(app)

with open('../data/gene_list.txt', 'r') as fh:
    entries = [l.strip().split(' ') for l in fh.readlines()]
    channel_labels = [(e[0], ('%s (%s)' % (e[0], e[1])
                              if len(e) > 1 else e[0]))
                      for e in entries]


def sorted_json_string(json_thing):
    """Produce a string that is unique to a json's contents."""
    if isinstance(json_thing, str):
        return json_thing
    elif isinstance(json_thing, (tuple, list)):
        return '[%s]' % (','.join(sorted(sorted_json_string(s)
                                         for s in json_thing)))
    elif isinstance(json_thing, dict):
        return '{%s}' % (','.join(sorted(k + sorted_json_string(v)
                                         for k, v in json_thing.items())))
    elif isinstance(json_thing, (int, float)):
        return str(json_thing)
    else:
        raise TypeError('Invalid type: %s' % type(json_thing))


def _get_query_hash(query_json):
    """Create an FNV-1a 32-bit hash from the query json"""
    return fnv1a_32(sorted_json_string(query_json).encode('utf-8'))



query_cache = {}


class ChannelSearchForm(FlaskForm):
    inhibits = SelectMultipleField(label='inhibit all of...',
                                   id='inhibit-select',
                                   choices=channel_labels)
    not_inhibits = SelectMultipleField(label='and not inhibit any of...',
                                       id='not-inhibit-select',
                                       choices=channel_labels)
    submit_button = SubmitField('Search')


@app.route('/get_reg_stmts', methods=['POST'])
def get_reg_stmts():
    pass


@app.route('/channel_search', methods=['POST'])
def channel_search():
    query = {}

    # Get inhibits
    inhibits = request.form.getlist('inhibits')
    query['inhibits'] = inhibits

    # Get not inhibits
    not_inhibits = request.form.getlist('not_inhibits')
    query['not_inhibits'] = not_inhibits

    agents = find_reg_non_regs(inhibits, not_inhibits)
    regs = []
    for agent in agents:
        urls = get_agent_urls(agent)
        regs.append((agent.name, urls))

    # Get results
    query['response_list'] = regs

    query_hash = _get_query_hash(query)
    if query_hash not in query_cache:
        query_cache[query_hash] = query
    session['query_hash'] = query_hash
    return redirect('/')


@app.route('/')
def index():
    channel_search_form = ChannelSearchForm()
    kwargs = {'channel_search_form': channel_search_form}
    if session.get('query_hash'):
        query_hash = session.pop('query_hash')
        response_list = query_cache[query_hash].get('response_list', [])
        old_inhibits_list = query_cache[query_hash].get('inhibits', [])
        old_not_inhibits_list = query_cache[query_hash].get('not_inhibits', [])
    else:
        response_list = []
        old_inhibits_list = []
        old_not_inhibits_list = []
    kwargs['response_list'] = response_list
    kwargs['old_inhibits_list'] = old_inhibits_list
    kwargs['old_not_inhibits_list'] = old_not_inhibits_list
    return render_template('index.html', **kwargs)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        port = int(sys.argv[1])
    else:
        port = 5000
    app.run(host='0.0.0.0', port=port)
