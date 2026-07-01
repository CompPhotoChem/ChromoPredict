#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Packaging for the ChromoPredict toolbox.

Installing is *not* strictly required -- the scripts also run straight from a
checkout if you put ``src/`` on the path (see ``-chromopredict-src`` in
``chromopredict_batch.py``).  But an editable install makes ``import
chromopredict`` work everywhere (tutorial notebooks, the batch script, ...)::

    pip install -e .

    python chromopredict_batch.py -input structures.xlsx -out predictions.csv \
           -plot error.png -rules auto woodward woodward_refine -solvent CCO
"""

from setuptools import find_namespace_packages, setup

setup(
    name="chromopredict",
    version="1.1.0",
    description="Automated Woodward-Fieser / Fieser / Fieser-Kuhn lambda_max prediction",
    # src/ layout: tell setuptools the packages live under src/
    package_dir={"": "src"},
    packages=find_namespace_packages(where="src", include=["chromopredict*"]),
    python_requires=">=3.9",
    install_requires=["numpy", "pandas", "pillow"],
)
