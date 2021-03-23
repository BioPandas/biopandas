#!/usr/bin/env python

from setuptools import setup
import re
import os
import sys


long_description = (
    "MkDocs is a fast, simple and downright gorgeous static site generator "
    "that's geared towards building project documentation. Documentation "
    "source files are written in Markdown, and configured with a single YAML "
    "configuration file."
)


def get_version(package):
    """Return package version as listed in `__version__` in `init.py`."""
    init_py = open(os.path.join(package, '__init__.py')).read()
    return re.search("__version__ = ['\"]([^'\"]+)['\"]", init_py).group(1)


def get_packages(package):
    """Return root package and all sub-packages."""
    return [dirpath
            for dirpath, dirnames, filenames in os.walk(package)
            if os.path.exists(os.path.join(dirpath, '__init__.py'))]


if sys.argv[-1] == 'publish':
    if os.system("pip freeze | grep wheel"):
        print("wheel not installed.\nUse `pip install wheel`.\nExiting.")
        sys.exit()
    if os.system("pip freeze | grep twine"):
        print("twine not installed.\nUse `pip install twine`.\nExiting.")
        sys.exit()
    os.system("python setup.py sdist bdist_wheel")
    os.system("twine upload dist/*")
    print("You probably want to also tag the version now:")
    print("  git tag -a {0} -m 'version {0}'".format(get_version("mkdocs")))
    print("  git push --tags")
    sys.exit()


setup(
    name="mkdocs",
    version=get_version("mkdocs"),
    url='https://www.mkdocs.org',
    license='BSD',
    description='Project documentation with Markdown.',
    long_description=long_description,
    author='Tom Christie',
    author_email='tom@tomchristie.com',  # SEE NOTE BELOW (*)
    packages=get_packages("mkdocs"),
    include_package_data=True,
    install_requires=[
        'click>=3.3',
        'Jinja2>=2.10.1',
        'livereload>=2.5.1',
        'lunr[languages]==0.5.8',  # must support lunr.js version included in search
        'Markdown>=3.2.1',
        'PyYAML>=3.10',
        'tornado>=5.0',
        'ghp-import>=1.0',
        'pyyaml_env_tag>=0.1'
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'mkdocs = mkdocs.__main__:cli',
        ],
        'mkdocs.themes': [
            'mkdocs = mkdocs.themes.mkdocs',
            'readthedocs = mkdocs.themes.readthedocs',
        ],
        'mkdocs.plugins': [
            'search = mkdocs.contrib.search:SearchPlugin',
        ],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        'Topic :: Documentation',
        'Topic :: Text Processing',
    ],
    zip_safe=False,
)

# (*) Please direct queries to the discussion group:
#     https://groups.google.com/forum/#!forum/mkdocs
