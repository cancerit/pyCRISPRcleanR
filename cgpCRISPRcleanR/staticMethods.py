import sys
import os
import tarfile
import tempfile
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE, STDOUT
import shlex
import re
import time
import logging.config
#os.environ["LD_LIBRARY_PATH"] = "/software/R-3.3.0/lib/R/lib"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/R-3.3.0/lib/R/li

from . import segmentation

from beautifultable import BeautifulTable
pd.options.display.width = 200

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
json_config = configdir + 'fileTypes.json'
logging.config.fileConfig(log_config)

log = logging.getLogger('cgpCRISPRcleanR')


class StaticMthods(object):
    """ Static methosds for common tasks """

    def __init__(self):
        super().__init__()

    def input_checker(infile):
        """
          checks user input file and returns it's type
        """
        try:
            if tarfile.is_tarfile(infile):
                log.info(("input is an archive:", infile))
                return 'tar'
            else:
                log.info(("input is a file:", infile))
                return 'file'
        except IsADirectoryError:
            return 'dir'
        except IOError as ioe:
            sys.exit('Error in reading input file:{}'.format(ioe.args[0]))

    def create_dict(inputlist):
        chrdict={}
        chr_count=1
        for chr in inputlist:
            chrdict[chr]=chr_count
            chr_count+=1
        return chrdict

    def format_counts(countfile,libfile,controls,min_read_count=30):
        col_to_drop=['gene','EXONE','CODE','STRAND']
        rename_col_dict={'CHRM': 'CHR', 'STARTpos':'startp','ENDpos':'endp','GENES':'genes'}
        print("Creating dataframe from input file:".format(countfile))
        counts=pd.read_csv(countfile, sep="\t", index_col='sgRNA')
        libdata=pd.read_csv(libfile, sep="\t", index_col='sgRNA')
        cldf = pd.concat([counts,libdata], axis=1, join='inner') # union of of indexes to create combined df count and library data
       # perform some cleanup here [ change row names and drop duplicate columns]
        cldf.drop(col_to_drop, axis=1, inplace=True)
        cldf=cldf.rename(columns=rename_col_dict)
        cldf.drop(cldf[cldf.iloc[:,0:controls].mean(axis=1) < min_read_count].index, inplace=True)
        # done with cleanup do some real work .....
        # create normalized count
        normed=cldf.iloc[:,0:cldf.columns.get_loc('genes')].div(cldf.iloc[:,0:cldf.columns.get_loc('genes')].agg('sum'))*10e6
        # add normalised count to origina dataframe
        #cldf.iloc[:,0:cldf.columns.get_loc('genes')]=normed # use if raw count is notrequired, replaces existing counts columns
        cldf=cldf.join(normed,rsuffix='_nc')
        fc = normed.apply(lambda x: np.log2( (x+0.5)/(normed.iloc[:,0:controls].mean(axis=1)+0.5) ) )
        fc.drop(fc.columns[0:controls], axis=1, inplace=True)
        no_rep=len(fc.columns)
        cldf['avgFC'] = fc.mean(axis=1) # claculate mean folchage and add to main data frame containg libraty annoatations
        cldf['BP'] = round( cldf['startp'] +  (cldf['endp'] - cldf['startp'] ) / 2  ).astype(int)
        cldf.sort_values(by=['CHR','startp'], ascending=True, inplace=True)
        return (cldf,no_rep)

    def genomwide_clean_chr(**kwargs):
        fc=kwargs['logfc']
        min_genes=kwargs.get('minTargetedGenes', 3)
        correctedFC=pd.DataFrame()
        segments=pd.DataFrame()
        for chridx in fc.CHR.unique():
            if chridx in ("1",'2'):
                chrdf=fc[fc.CHR == chridx ]
                print(chridx)
                (corrected_fc,regions)=segmentation.do_segmentation(chrdf)
                correctedFC=correctedFC.append(corrected_fc)
                segments=segments.append(regions)
        # reindex data frame after concatenation
        segments.reset_index(drop=True, inplace=True)
        return (correctedFC,segments)

    def corrected_counts(alldf,segments,minTargetedGenes=3,no_rep=1,sample_id='HT-29',controls=1):
        #sgRNA <controls> <Sample Raw count>   genes CHR  startp    endp  <normalised count>   avgFC      BP  correction  correctedFC
        # store normalised counts,selected segmen locations will be overwiteen as reverted counts
        revertedCounts=alldf.iloc[:,alldf.columns.get_loc('endp')+controls+1:alldf.columns.get_loc('avgFC')]
        for seg in segments.itertuples():
            i=seg.Index
            start=seg.startRow
            end=seg.endRow
            ntarg=seg.nGenes
            idxs=list(range(start, end))
            if ntarg>=minTargetedGenes:
                #get a segment slice of a dataframe
                reverted=pd.DataFrame()
                segdata=alldf.iloc[idxs]
                #average fold change
                FC=segdata.avgFC
                #mean normalised control count
                nc=segdata.iloc[:,segdata.columns.get_loc('endp')+1:segdata.columns.get_loc('endp')+controls+1]
                c=nc.mean(axis=1)
                #print(segdata.head())
                N=segdata.correctedFC
                CF = N - FC
                reverted['revc'] = c*(pow(2,N))
                normed_num=segdata.iloc[:,segdata.columns.get_loc('endp')+controls+1:segdata.columns.get_loc('avgFC')]
                normed_num+=1
                proportions=normed_num.div(normed_num.agg('sum', axis=1),axis=0)
                reverted=reverted*no_rep
                revertedCounts.iloc[idxs]=proportions.mul(reverted.revc,axis=0)
        # store reverted count to original data frame
        alldf=alldf.join(revertedCounts,rsuffix='_rev')
        print(alldf.shape)
        return alldf




























































    @staticmethod
    def get_position(row):
        return round(row['startp'] + (row['endp'] - row['startp']) / 2)
    @staticmethod
    def swap_key2val(inputdict):
        return {val:key for key,val in inputdict.items()}
