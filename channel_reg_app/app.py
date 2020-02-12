from flask_wtf import Form
from flask_bootstrap import Bootstrap
from flask import Flask, render_template, request
from wtforms import SubmitField, StringField, validators, SelectMultipleField

app = Flask(__name__)
app.config['SECRET_KEY'] = 'dev_key'
Bootstrap(app)


channels = ['A', 'B', 'C']


class ChannelSearchForm(Form):
    inhibits = SelectMultipleField(label='Inhibits',
                                   choices=[(c, c) for c in channels],
                                   validators=[validators.unicode_literals])
    not_inhibits = SelectMultipleField(label='Does not inhibit',
                                       choices=[(c, c) for c in channels],
                                       validators=[validators.unicode_literals])
    submit_button = SubmitField('Search')


@app.route('/channel_search', methods=['POST'])
def channel_search():
    inhibits = request.form.getlist('inhibits')
    print(inhibits)
    not_inhibits = request.form.getlist('not_inhibits')
    print(not_inhibits)
    return str(inhibits) + str(not_inhibits)


@app.route('/')
def index():
    channel_search_form = ChannelSearchForm()
    kwargs = {'channel_search_form': channel_search_form}
    return render_template('index.html', **kwargs)


if __name__ == '__main__':
    app.run()
