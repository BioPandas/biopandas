from os.path import realpath, dirname, join
from setuptools import setup, find_packages
import biopandas

VERSION = biopandas.__version__
PROJECT_ROOT = dirname(realpath(__file__))

REQUIREMENTS_FILE = join(PROJECT_ROOT, 'requirements.txt')

with open(REQUIREMENTS_FILE) as f:
    install_reqs = f.read().splitlines()

install_reqs.append('setuptools')

setup(name='biopandas',
      version=VERSION,
      description='Machine Learning Library Extensions',
      author='Sebastian Raschka',
      author_email='mail@sebastianraschka.com',
      url='https://github.com/rasbt/biopandas',
      packages=find_packages(),
      package_data={'': ['LICENSE.txt',
                         'README.md',
                         'requirements.txt']
                    },
      include_package_data=True,
      install_requires=install_reqs,
      license='BSD 3-Clause',
      platforms='any',
      classifiers=[
             'License :: OSI Approved :: BSD License',
             'Development Status :: 5 - Production/Stable',
             'Operating System :: Microsoft :: Windows',
             'Operating System :: POSIX',
             'Operating System :: Unix',
             'Operating System :: MacOS',
             'Programming Language :: Python :: 3',
             'Programming Language :: Python :: 3.5',
             'Programming Language :: Python :: 3.6',
             'Programming Language :: Python :: 3.7',
             'Topic :: Scientific/Engineering',
      ],
      long_description="""
Biopandas is a python package for working with molecular structures
in pandas DataFrames.


Contact
=============

If you have any questions or comments about biopandas,
please feel free to contact me via
eMail: mail@sebastianraschka.com
or Twitter: https://twitter.com/rasbt

This project is hosted at https://github.com/rasbt/biopandas

The documentation can be found at http://rasbt.github.io/biopandas/

""")
