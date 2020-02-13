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


class ChannelSearchForm(FlaskForm):
    inhibits = SelectMultipleField(label='Inhibits',
                                   id='inhibit-select',
                                   choices=channel_labels,
                                   validators=[validators.unicode_literals])
    not_inhibits = SelectMultipleField(label='Does not inhibit',
                                       id='not-inhibit-select',
                                       choices=channel_labels,
                                       validators=[validators.unicode_literals])
    submit_button = SubmitField('Search')


@app.route('/channel_search', methods=['POST'])
def channel_search():
    # Get inhibits
    inhibits = request.form.getlist('inhibits')
    session['inhibits'] = inhibits

    # Get not inhibits
    not_inhibits = request.form.getlist('not_inhibits')
    session['not_inhibits'] = not_inhibits

    agents = find_reg_non_regs(inhibits, not_inhibits)
    regs = []
    for agent in agents:
        urls = get_agent_urls(agent)
        agent_str = '%s %s' % (agent.name, ', '.join(urls))
        regs.append(agent_str)

    # Get results
    session['response_list'] = regs
    return redirect('/')


@app.route('/')
def index():
    channel_search_form = ChannelSearchForm()
    kwargs = {'channel_search_form': channel_search_form}
    if session.get('response_list'):
        response_list = session.pop('response_list', [])
        old_inhibits_list = session.pop('inhibits', [])
        old_not_inhibits_list = session.pop('not_inhibits', [])
    else:
        response_list = []
        old_inhibits_list = []
        old_not_inhibits_list = []
    kwargs['response_list'] = response_list
    kwargs['old_inhibits_list'] = old_inhibits_list
    kwargs['old_not_inhibits_list'] = old_not_inhibits_list
    return render_template('index.html', **kwargs)


if __name__ == '__main__':
    app.run()
