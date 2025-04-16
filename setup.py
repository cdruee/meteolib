#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

from setuptools import setup

v = {}
v_path = os.path.join(*'meteolib/_version.py'.split('/'))
with open(v_path) as v_file:
    exec(v_file.read(), v)

setup(
    name=v['__title__'],
    version=v['__version__'],
    packages=[v['__title__']],
    package_dir={v['__title__']: v['__title__']},
    license=v['__license__'],
    author=v['__author__'],
    author_email=v['__author_email__'],
    url=v['__url__'],
    description=v['__description__'],
    long_description=open('README.md').read(),
    install_requires=[
        'numpy',
        'pandas',
    ],
    keywords='meteorology',
)
