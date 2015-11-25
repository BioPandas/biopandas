from setuptools import setup, find_packages

def calculate_version():
    initpy = open('biopandas/__init__.py').read().split('\n')
    version = list(filter(lambda x: '__version__' in x, initpy))[0].split('\'')[1]
    return version

package_version = calculate_version()

setup(name='biopandas',
      version=package_version,
      description='Molecular Structures in Pandas DataFrames',
      author='Sebastian Raschka',
      author_email='mail@sebastianraschka.com',
      url='https://github.com/rasbt/biopandas',
      license='new BSD',
      zip_safe=True,
      packages=find_packages(),
      platforms='any',
      install_requires=['numpy', 'pandas'],
      keywords=['bioinformatics', 'molecular structures',
                'protein databank', 'computational biology',
                'protein structures'],
      classifiers=[
             'License :: OSI Approved :: BSD License',
             'Development Status :: 5 - Production/Stable',
             'Operating System :: Microsoft :: Windows',
             'Operating System :: POSIX',
             'Operating System :: Unix',
             'Operating System :: MacOS',
             'Programming Language :: Python :: 2',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Python :: 3',
             'Programming Language :: Python :: 3.3',
             'Programming Language :: Python :: 3.4',
             'Programming Language :: Python :: 3.5',
             'Topic :: Scientific/Engineering',
      ],
      long_description="""
Biopandas is a python package for working with molecular structures
in pandas DataFrames.


Contact
=============

If you have any questions or comments about biopandas, please feel free to contact me via
eMail: mail@sebastianraschka.com
or Twitter: https://twitter.com/rasbt

This project is hosted at https://github.com/rasbt/biopandas

The documentation can be found at http://rasbt.github.io/biopandas/

""",
    )
