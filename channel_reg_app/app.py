from flask_wtf import Form
from flask_bootstrap import Bootstrap
from flask import Flask, render_template, request, redirect, session
from wtforms import SubmitField, StringField, validators, SelectMultipleField

app = Flask(__name__)
app.config['SECRET_KEY'] = 'dev_key'

Bootstrap(app)


channels = ['A', 'B', 'C']


class ChannelSearchForm(Form):
    inhibits = SelectMultipleField(label='Inhibits',
                                   id='inhibit-select',
                                   choices=[(c, c) for c in channels],
                                   validators=[validators.unicode_literals])
    not_inhibits = SelectMultipleField(label='Does not inhibit',
                                       id='not-inhibit-select',
                                       choices=[(c, c) for c in channels],
                                       validators=[validators.unicode_literals])
    submit_button = SubmitField('Search')


@app.route('/channel_search', methods=['POST'])
def channel_search():
    # Get inhibits
    inhibits = request.form.getlist('inhibits')
    print(inhibits)
    session['inhibits'] = inhibits

    # Get not inhibits
    not_inhibits = request.form.getlist('not_inhibits')
    session['not_inhibits'] = not_inhibits
    print(not_inhibits)

    # Get results
    response_list = list(set(inhibits).union(set(not_inhibits)))
    session['response_list'] = response_list
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
