import io
import os
from setuptools import setup, find_packages
from setuptools.command.install import install
import urllib.request
import shutil
import sys
import zipfile
import tarfile

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
        print('Installing TMalign compiled code')
        target_folder = os.path.join(PROJECT_ROOT, 'biopandas', 'align')
        os.makedirs(target_folder, exist_ok=True)
        
        # Define the file names for each OS with format info
        files = {
            'win32': ('USalignWin64.zip', 'zip'),
            'linux': ('USalignLinux64.zip', 'zip'),
            'darwin': ('USalignMac.tar.gz', 'tar.gz')  # Fixed: Mac uses tar.gz
        }

        # Detect the OS
        os_name = sys.platform
        print(f'DEBUG: Detected OS: {os_name}')
        print(f'DEBUG: Available files: {files}')

        # Select the appropriate file based on the OS
        if os_name in files:
            selected_file, file_format = files[os_name]
            print(f'DEBUG: Selected file: "{selected_file}" (format: {file_format})')
        else:
            print(f'ERROR: Unsupported OS: {os_name}')
            raise Exception('Unsupported OS')

        # Download the file
        url = f'https://zhanggroup.org/US-align/bin/module/{selected_file}'
        print(f'DEBUG: Final URL: {url}')
        
        try:
            print(f'DEBUG: Starting download...')
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'})
            with urllib.request.urlopen(req) as response, open(selected_file, 'wb') as out_file:
                out_file.write(response.read())
            print(f'DEBUG: Download completed successfully')
            print(f'DEBUG: File size: {os.path.getsize(selected_file)} bytes')
        except Exception as e:
            print(f'ERROR: Download failed: {e}')
            print(f'ERROR: URL was: {url}')
            raise

        # Extract file based on format
        if file_format == 'zip':
            extract_dir = selected_file.replace(".zip", "")
            print(f'DEBUG: Extracting ZIP to: {extract_dir}')
            with zipfile.ZipFile(selected_file, "r") as zip_ref:
                zip_ref.extractall(extract_dir)
        else:  # tar.gz
            extract_dir = selected_file.replace(".tar.gz", "")
            print(f'DEBUG: Extracting TAR.GZ to: {extract_dir}')
            with tarfile.open(selected_file, "r:gz") as tar_ref:
                tar_ref.extractall(extract_dir)

        # Move the compiled code to the target folder
        print(f'DEBUG: PROJECT_ROOT: {PROJECT_ROOT}')
        print(f'DEBUG: os_name: {os_name}')
        print(f'DEBUG: extract_dir: {extract_dir}')
        print(f'DEBUG: target_folder: {target_folder}')
        
        unzipped_path = os.path.join(extract_dir, 'USalign')
        print(f'DEBUG: unzipped_path: {unzipped_path}')
        
        if os_name == 'win32':
            source_file = os.path.join(unzipped_path, 'USalign.exe')
            target_file = os.path.join(target_folder, 'USalign.exe')
        else:
            source_file = os.path.join(unzipped_path, 'USalign')
            target_file = os.path.join(target_folder, 'USalign')
        
        print(f'DEBUG: source_file: {source_file}')
        print(f'DEBUG: target_file: {target_file}')
        print(f'DEBUG: source_file exists: {os.path.exists(source_file)}')
        
        if os.path.exists(source_file):
            os.replace(source_file, target_file)
            print(f'SUCCESS: TMalign installed at {target_file}')
        else:
            print(f'ERROR: Source file not found: {source_file}')
            # List directory contents for debugging
            if os.path.exists(extract_dir):
                print(f'DEBUG: Contents of {extract_dir}:')
                for root, dirs, files in os.walk(extract_dir):
                    level = root.replace(extract_dir, '').count(os.sep)
                    indent = ' ' * 2 * level
                    print(f'{indent}{os.path.basename(root)}/')
                    subindent = ' ' * 2 * (level + 1)
                    for f in files:
                        print(f'{subindent}{f}')

        # Remove the downloaded and extracted files
        shutil.rmtree(extract_dir)
        os.remove(selected_file)
        print('DEBUG: Cleanup completed')

setup(name='biopandas',
      version=VERSION,
      description='Machine Learning Library Extensions',
      author='Sebastian Raschka',
      author_email='mail@sebastianraschka.com',
      url='https://github.com/rasbt/biopandas',
      packages=find_packages(),
      package_data={'': ['LICENSE.txt',
                         'README.md',
                         'requirements.txt',
                         'USalign.exe',
                         'USalign']
                    },
      include_package_data=True,
      install_requires=install_reqs,
      extras_require={'test': ['pytest', 'pytest-cov','flake8'],},
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
      cmdclass={'install': TmAlignInstall},
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
