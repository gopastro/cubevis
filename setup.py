#!/usr/bin/env python
from setuptools import setup, find_packages

NAME = 'cubevis'
VERSION = '0.1'

setup(
    name=NAME,
    version=VERSION,
    description='Spectral Line Data Cube Reduction and Visualization Tool for Radio Astronomy',
    author='Gopal Narayanan <gopal@astro.umass.edu>',
    packages=find_packages(),
    )
