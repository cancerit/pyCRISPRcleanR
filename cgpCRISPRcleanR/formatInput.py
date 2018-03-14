import csv
import numpy as np
import pandas as pd
import os
import sys
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
        """ check input type and presence of user supplied
          input files
        """
        super().check_input()
        input_type = []
        for infile in (self.countfile, self.libfile):
            input_type.append(sm.input_checker(infile))
            print("File:{} is type:{}".format(infile,input_type))
        return None
    def run_analysis(self):
        """
          method to run the analysis
        """
        controls=self.ncontrols
        min_read_count=self.minreads
        min_target_genes=self.mingenes
        ignoredGenes=[] # should come from command line ??? 
        sample='mysample'
        cpus=10
        #check input files
        self.check_input()
        cldf = sm.combine_count_n_library(self.countfile,self.libfile)
        cldf = sm.cleanup_data(cldf)
        cldf = sm.filter_data(cldf,controls,min_read_count)
        (cldf,num_rep)=sm.get_norm_count_n_fold_changes(cldf,controls)
        cbs_dict = sm.run_cbs(cldf,cpus,sample)
        all_data = sm.process_segments(cbs_dict,ignoredGenes,min_target_genes,controls,num_rep)
        all_data.to_csv("corrected_counts_alldf_v9nc.tsv", sep='\t', mode='w', header=True, index=True, index_label='gRNA',doublequote=False)
        print(all_data[all_data.gene=='A1BG'])
