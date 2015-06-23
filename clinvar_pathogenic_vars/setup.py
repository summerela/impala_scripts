#!/usr/bin/env python

#installation file
from setuptools import setup

setup(name='clinvar_pathogenic_vars',
      version='0.1',
      description='Find potentially pathogenic variants using impala and clinvar',
      platforms=['Windows', 'Linux', 'OSX'],
      url='https://github.com/summerela/impala_scripts/clinvar_pathogenic_vars',
      author='Summer Elasady',
      author_email='selasady@systemsbiology.org',
      license='ISB',
      packages=['clinvar_query'],
      install_requires=['setuptools'],
      dependency_links=['http://pypi.python.org/pypi/pandas/0.16.2/#downloads',
                        'https://pypi.python.org/packages/source/i/impyla/impyla-0.10.0.tar.gz',
                        'https://pypi.python.org/packages/source/s/sqlparse/sqlparse-0.1.15.tar.gz'
                        ],
      include_package_data=True,
      zip_safe=False,
      scripts=['clinvar_query.py'])