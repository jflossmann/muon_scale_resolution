import json
import yaml
import os
import sys

dname = sys.argv[1]
datasets = f'data/{dname}'

# open file with dataset names
with open(datasets) as f:
    dsets = yaml.load(f, yaml.Loader)

lib = "root://xrootd-cms.infn.it//"
fdict = {}

for k in dsets.keys():
    for d in dsets[k]["names"]:
        print(f'dasgoclient -query="file dataset={d}"')
        stream = os.popen(f'dasgoclient -query="file dataset={d}"')
        fdict[k] = [
            lib+s.replace('\n', '') for s in stream.readlines()
        ]

with open(datasets.replace("datasets", 'nanoAODs'), "w") as f:
    yaml.dump(fdict, f)