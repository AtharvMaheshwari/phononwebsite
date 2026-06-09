import json
import gzip
from python.phononweb import Phonon
from python.phononweb.chirality import compute_chiral_properties

with gzip.open('data/localdb/mp-3427.json.gz', 'rt') as f:
    data = json.load(f)

ph = Phonon()
ph.read_json(data)

# This will compute everything
chiral_data = ph.compute_chiral_properties()

# Inject the calculated properties
data.update(chiral_data)

# Write to the new Al4Li4O8 directory
with open('data/localdb/al4li4o8/data.json', 'w') as f:
    json.dump(data, f)
with open('build/data/localdb/al4li4o8/data.json', 'w') as f:
    json.dump(data, f)

# print("Actual properties calculated and saved!")