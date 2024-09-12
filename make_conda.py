import yaml
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '')] + sys.path

import src.version as version

data = {
    "package": {
        "name": "cgphylo",
        "version": version.__version__
    },
    "source": {
        "url": "https://github.com/genomicepidemiology/CGPhylo/archive/refs/tags/{}.tar.gz".format(version.__version__),
    },
    "build": {
        "number": 0,
        "noarch": "python",
        "script": "{{ PYTHON }} -m pip install . --no-deps --ignore-installed -vvv",
        "script_env": ["PYTHONNOUSERSITE=1"]
    },
    "requirements": {
        "host": [
            "python >=3.6",
            "pip"
        ],
        "run": [
            "kma >=1.4.9",
            "biopython"
        ]
    },
    "about": {
        "home": "https://github.com/genomicepidemiology/CGPhylo",
        "summary": "CGPhylo - core genome phylogenetics",
        "license": "Apache-2.0"
    },
    "extra": {
        "recipe-maintainers": ["malhal"]
    }
}

# Convert the data to YAML and print it
os.system('mkdir -p conda')
yaml_str = yaml.dump(data, sort_keys=False)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)
