#!/usr/bin/env python3

from setuptools import setup

config = {
    'version': '1.0.0',
    'name': 'cgpCRISPRcleanR',
    'description': 'tool to comapre files and/or archives',
    'author': 'Shriram G Bhosle',
    'url': 'https://github.com/CancerIT/cgpCRISPRcleanR',
    'author_email': 'cgphelp@sanger.ac.uk',
    'python_requires': '>= 3.3',
    'setup_requires': ['pytest','pytest-cov'],
    'install_requires': ['logging','rpy2','pandas'],
    'packages': ['cgpCRISPRcleanR'],
    'package_data': {'cgpCRISPRcleanR': ['config/*.json','config/*.conf']},
    'entry_points': {
        'console_scripts': ['cgp_crispr_cleanr=cgpCRISPRcleanR.cgpcrisprCleanR_command:main'],
    }
}

setup(**config)
