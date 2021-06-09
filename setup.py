#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pysodist",
    version="0.0.5",
    author="Joey Davis",
    author_email="jhdavis@mit.edu",
    description="A python-based implementation of the isodist fitting routines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jhdavislab/pysodist",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'pysodist = pysodist.__main__:main',
        ],
    },
    include_package_data=True,
    python_requires='>=3.6',
    install_requires=[
        'matplotlib>=3.3.2',
        'scipy>=1.5.2',
        'seaborn>=0.11.0',
        'pandas>=1.1.3',
        'numpy>=1.19.2',
        'pyteomics>=4.3.3',
        'qgrid>=1.3.1']

)