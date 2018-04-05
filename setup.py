#!/usr/bin/env python3

from setuptools import setup

config = {
    'version': '1.0.0',
    'name': 'pyCRISPRcleanR',
    'description': 'tool to comapre files and/or archives',
    'author': 'Shriram G Bhosle',
    'url': 'https://github.com/CancerIT/pyCRISPRcleanR',
    'author_email': 'pyhelp@sanger.ac.uk',
    'python_requires': '>= 3.3',
    'setup_requires': ['pytest','pytest-cover'],
    'install_requires': ['rpy2','pandas', 'numpy'],
    'packages': ['pyCRISPRcleanR'],
    'package_data': {'pyCRISPRcleanR':['config/*.conf']},
    'entry_points': {
        'console_scripts': ['pyCRISPRCleanR=pyCRISPRcleanR.crisprCleanR_command:main'],
    }
}

setup(**config)
