#!/usr/bin/env python3

from setuptools import setup

config = {
    'version': '2.0.7',
    'name': 'pyCRISPRcleanR',
    'description': 'This is python implementation of CRISPRcleanR package for unsupervised identification and correction of gene independent cell responses to CRISPR-cas9 targeting',
    'author': 'Shriram Bhosle',
    'url': 'https://github.com/CancerIT/pyCRISPRcleanR',
    'author_email': 'cgphelp@sanger.ac.uk',
    'python_requires': '>= 3.3',
    'setup_requires': ['pytest','pytest-cover', 'radon'],
    'install_requires': ['scipy','rpy2', 'pandas', 'numpy', 'plotly', 'tzlocal'],
    'packages': ['pyCRISPRcleanR'],
    'package_data': {'pyCRISPRcleanR':['config/*.conf','config/*.json','config/*.tar.gz','segmentation/*.py']},
    'entry_points': {
        'console_scripts': ['pyCRISPRcleanR=pyCRISPRcleanR.crisprCleanR_command:main'],
    }
}

setup(**config)
