import io
import os
from setuptools import setup, find_packages
from setuptools.command.install import install
import urllib.request
import shutil
import sys
import zipfile

VERSION = None
with io.open(
    os.path.join(os.path.dirname(__file__), 'biopandas/__init__.py'),
    encoding='utf-8'
) as f:
    for l in f:
        if not l.startswith('__version__'):
            continue
        VERSION = l.split('=')[1].strip(' "\'\n')
        break
PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__))

REQUIREMENTS_FILE = os.path.join(PROJECT_ROOT, 'requirements.txt')

with open(REQUIREMENTS_FILE) as f:
    install_reqs = f.read().splitlines()

install_reqs.append('setuptools')

class TmAlignInstall(install):
    def run(self):
        install.run(self)
        target_folder = os.path.join(PROJECT_ROOT, 'biopandas', 'align')
        os.makedirs(target_folder, exist_ok=True)  # Ensure the subfolder exists
        # download the TMalign compiled code
        # Define the file names for each OS
        files = {
            'win32': 'USalignWin64.zip',
            'linux': 'USalignLinux64.zip',
            'darwin': 'USalignMac.zip'
        }

        # Detect the OS and architecture
        os_name = sys.platform

        # Select the appropriate file based on the OS
        if os_name in files:
            selected_file = files[os_name]
        else:
            raise Exception('Unsupported OS')

        # Download the file
        url = f'https://zhanggroup.org/US-align/bin/module/{selected_file}'
        urllib.request.urlretrieve(url, selected_file)

        # unzip the file
        with zipfile.ZipFile(selected_file, "r") as zip_ref:
            zip_ref.extractall(selected_file.replace(".zip", ""))

        # Mv the compiled code to the target folder. It is either USalign or USalign.exe on windows
        target_folder = os.path.join(PROJECT_ROOT, 'biopandas', 'align')
        unzipped_path = f'{selected_file.replace(".zip", "")}/USalign/'
        if os_name == 'win32':
            os.replace(os.path.join(unzipped_path, 'USalign.exe'), os.path.join(target_folder, 'USalign.exe'))
        else:
            # on linux and mac
            os.replace(os.path.join(unzipped_path, 'USalign'), os.path.join(target_folder, 'USalign'))

        # Remove the downloaded and the unzipped files (everything starting with USalign)
        shutil.rmtree(selected_file.replace(".zip", ""))
        os.remove(selected_file)

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
      extras_require={'test': ['pytest', 'pytest-cov','flake8', 'nose'],},
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
             'Programming Language :: Python :: 3.8',
             'Programming Language :: Python :: 3.9',
             'Topic :: Scientific/Engineering',
      ],
      cmdclass={'install': TmAlignInstall}, # for handling TMalign compilation
      long_description_content_type='text/markdown',
      long_description="""
Biopandas is a Python package for working with molecular structures
in pandas DataFrames.


Contact
=============

If you have any questions or comments about biopandas,
please feel free to contact me via
eMail: mail@sebastianraschka.com
or Twitter: https://twitter.com/rasbt

This project is hosted at https://github.com/rasbt/biopandas

The documentation can be found at http://rasbt.github.io/biopandas/

"""
)
