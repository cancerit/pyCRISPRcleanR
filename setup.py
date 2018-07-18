#!/usr/bin/env python3

from setuptools import setup

config = {
    'version': '1.1.2',
    'name': 'pyCRISPRcleanR',
    'description': 'This is python implementation of CRISPRcleanR package for unsupervised identification and correction of gene independent cell responses to CRISPR-cas9 targeting',
    'author': 'Shriram Bhosle',
    'url': 'https://github.com/CancerIT/pyCRISPRcleanR',
    'author_email': 'cgphelp@sanger.ac.uk',
    'python_requires': '>= 3.3',
    'setup_requires': ['pytest','pytest-cover', 'radon'],
    'install_requires': ['rpy2', 'pandas', 'numpy', 'plotly', 'tzlocal'],
    'find_links' : 'https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.7.tar.gz/download',
    'packages': ['pyCRISPRcleanR'],
    'package_data': {'pyCRISPRcleanR':['config/*.conf','segmentation/*.py']},
    'entry_points': {
        'console_scripts': ['pyCRISPRCleanR=pyCRISPRcleanR.crisprCleanR_command:main'],
    }
}

setup(**config)
