import pickle
import sys

ALL_ENZYMES = sys.argv[1:]


print(ALL_ENZYMES)

with open(ALL_ENZYMES[0], 'rb') as fh:
	en1 = pickle.load(fh)


with open(ALL_ENZYMES[1], 'rb') as fh:
	en2 = pickle.load(fh)

print(en1)
print(en2)
print(ALL_ENZYMES)
"""
OUTPUT1 = "/Users/sbunga/PycharmProjects/nextflow/output/pkl/run"
OUTPUT2 = "test"
dump_out = 'hello world'
with open(ALL_ENZYMES, 'rb') as fh:
	all_enzymes = pickle.load(fh)


with open(OUTPUT1, 'wb') as fh:
	pickle.dump(dump_out, fh)

with open(OUTPUT2, 'wb') as fh:
	pickle.dump(dump_out, fh)
"""