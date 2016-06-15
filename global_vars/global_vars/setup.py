#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='global_vars',
      version='1.0',
      description='Create table of all variants found, annotate and run snpeff for functional annot',
      long_description=readme(),
      url='https://github.com/summerela/impala_scripts/global_vars',
      author='Summer Rae Elasady',
      author_email='selasady@systemsbiology.net',
      license='ISB',
      packages=['global_vars'],
      install_requires=[
          'pandas', 'impyla',
      ], include_package_data=True,
      zip_safe=False)