import csv
import numpy as np
import pandas as pd
import json
import tarfile
import os
import sys
import tempfile
import math
import re
import matplotlib
import logging.config
from sys import stderr
from cgpCRISPRcleanR.abstractCrispr import AbstractCrispr
from cgpCRISPRcleanR.staticMethods import StaticMthods as sm

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
json_config = configdir + 'fileTypes.json'
logging.config.fileConfig(log_config)

log = logging.getLogger('cgpCRISPRcleanR')

'''
Compare file, folder or archived data in .tar, .gz, .bz2 format
input any of the above
output -- list of items that are common in two compared archives
based on the MD5sum and/or data contents
'''

class CrisprCleanR(AbstractCrispr):
    """
      Main class implements abstract class and
      its methods to check the inputype of a given file and loads parameters from user config
      file required in json format
    """

    def check_input(self):
        """ check input type of user supplied
          input filea and fileb keyword arguments
          and returns tupule of input file type
        """
        super().check_input()
        input_type = []
        for infile in (self.countfile, self.libfile):
            input_type.append(sm.input_checker(infile))
            print("File:{} is type:{}".format(infile,input_type) )
            if(infile == 'HT-29_counts.tsv'):
                (lib_count_n_fc,no_rep)=sm.format_counts(self.countfile,self.libfile,self.ncontrols,min_read_count=30)
                (correctedFC,segments)=sm.genomwide_clean_chr(logfc=lib_count_n_fc,minTargetedGenes=3)
                all_data=sm.corrected_counts(correctedFC,segments,minTargetedGenes=3,no_rep=no_rep)

        return tuple(input_type)
