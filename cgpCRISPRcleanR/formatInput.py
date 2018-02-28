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
                (lib_count_n_fc)=sm.format_counts(self.countfile,self.libfile,self.ncontrols,min_read_count=30)
                #(logFC)=sm.sortbyposFC(fc,grnalib,chrdict)
                #(correctedFC,segments)=sm.genomwide_clean_chr(logfc=logFC,min_genes=3)
                #sm.corrected_counts(nc,correctedFC,segments,grnalib,min_genes=3)

        return tuple(input_type)

    def format_counts(countfile,libfile,controls,min_read_count=30):
        col_to_drop=['gene','EXONE','CODE','STRAND']
        rename_col_dict={'CHRM': 'CHR', 'STARTpos':'startp','ENDpos':'endp','GENES':'genes'}

        print("Creating dataframe from input file:".format(countfile))
        counts=pd.read_csv(countfile, sep="\t", index_col='sgRNA')
        libdata=pd.read_csv(libfile, sep="\t", index_col='sgRNA')
        print(counts.head())
        print(libdata.head())
        counts_n_libdata = pd.concat([counts,libdata], axis=1, join='inner') # union of of indexes
        counts_n_libdata.drop(col_to_drop, axis=1, inplace=True)
        counts_n_libdata=counts_n_libdata.rename(columns=rename_col_dict)
        print(counts_n_libdata.head())
        # perform some cleanup here [ change row names and drop duplicate columns]

        chrdict=sm.create_dict(libdata.CHRM.unique())
        # process data for matched sgRNA's in library file
        matched_count=counts.loc[libdata.index]
        # index matched count data as per libdata index
        #matched_count.reindex(libdata.index, axis='index')
        # get only counts data start from column 1 as index column is not counted ??
        numd = matched_count.iloc[:,1:]
        # filter datapoints based on min read counts in control sample ???
        #print(numd.iloc[0:8,:])
        mask=numd.iloc[:,0:controls].mean(axis=1) >= min_read_count # asix 1 refers along the column aggregation
        numd=numd.loc[mask]
        counts=counts.loc[mask] # only used for visualization in future
        libdata=libdata.loc[mask]
        # check if index isidentical
        cidx=counts.index
        nidx=numd.index
        libdex=libdata.index
        if not cidx.difference(libdex).tolist():
            print("indexes are  matching")
            normed=numd.div(numd.agg('sum'))*10e6
            #print(normed.head())
            # get foldchage data
            foldchanges = normed.apply(lambda x: np.log2( (x+0.5)/(normed.iloc[:,0:controls].mean(axis=1)+0.5) ) )
            # drop control columns
            print("foldchanges Before:{} ".format(foldchanges.shape))
            foldchanges.drop(foldchanges.columns[0:controls], axis=1, inplace=True)
            print("numd:{} counts:{} libdata:{}: foldchanges:{} ".format(numd.shape,counts.shape,libdata.shape,foldchanges.shape))
            #print(foldchanges.head())
            return (foldchanges,normed,libdata,chrdict)
            #c_foldchanges<-log2((normed[,3+i]+0.5)/(normed[,3]+0.5))
        else:
            print("Unmatched indexes found:{}".format(cidx.difference(nidx).tolist()) )
