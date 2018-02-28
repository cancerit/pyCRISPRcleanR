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
pd.options.display.width = 220

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
        print(cldf.shape)
        cldf.drop(cldf[cldf.iloc[:,0:controls].mean(axis=1) < min_read_count].index, inplace=True)
        print(cldf.shape)
        #print(cldf.head())
        # done with cleanup do some real work .....
        # create normalized count
        normed=cldf.iloc[:,0:cldf.columns.get_loc('genes')].div(cldf.iloc[:,0:cldf.columns.get_loc('genes')].agg('sum'))*10e6
        fc = normed.apply(lambda x: np.log2( (x+0.5)/(normed.iloc[:,0:controls].mean(axis=1)+0.5) ) )
        fc.drop(fc.columns[0:controls], axis=1, inplace=True)
        cldf['avgFC'] = fc.mean(axis=1) # claculate mean folchage and add to main data frame containg libraty annoatations
        cldf['BP'] = round( cldf['startp'] +  (cldf['endp'] - cldf['startp'] ) / 2  ).astype(int)
        cldf.sort_values(by=['CHR','startp'], ascending=True, inplace=True)
        # joining two data frames keeping indeticalcolumns with added prefix
        #cldf=cldf.join((cldf.iloc[:,0:cldf.columns.get_loc('genes')].div(cldf.iloc[:,0:cldf.columns.get_loc('genes')].agg('sum'))*10e6), rsuffix='_nc' )
        print(cldf.head())
        return (cldf)

    def sortbyposFC(fc,grnalib,chrdict):
        fc['avgFC']=fc.mean(axis=1) # claculate mean folchage
        # drop all clumns except avgFC
        fc.drop(fc.columns[0:-1], axis=1, inplace=True)
        #print(grnalib.head())
        # discuss data format or explain requirement in README
        #libcols=grnalib.columns.tolist()
        libcols=['CHRM','STARTpos','ENDpos','GENES']
        grnalib=grnalib[libcols]
        annot_df=pd.concat([grnalib,fc], axis=1)
        annot_df=annot_df.rename(columns={'CHRM': 'CHR', 'STARTpos':'startp','ENDpos':'endp','GENES':'genes'})
        #annot_df['BP'] = annot_df.apply(StaticMthods.get_position,axis=1) #  takes more time to run
        annot_df['BP']= round( annot_df['startp'] +  (annot_df['endp'] - annot_df['startp'] ) / 2  ).astype(int)
        annot_df=annot_df.sort_values(by=['CHR','startp'], ascending=True)
        # add chromsome index to data frame
        #annot_df['CHR'].replace(chrdict, inplace=True) # uncomment if want to use indexed chr names
        #index_to_chr_dict=StaticMthods.swap_key2val(chrdict)
        #annot_df['CHR'].replace(index_to_chr_dict, inplace=True) # decoding back to original chr names
        return annot_df

    def genomwide_clean_chr(**kwargs):
        fc=kwargs['logfc']
        min_genes=kwargs.get('min_genes', 3)
        correctedFC=pd.DataFrame()
        segments=pd.DataFrame()
        for chridx in fc.CHR.unique():
            if chridx in ("1",'2'):
                chrdf=fc[fc.CHR == chridx ]
                print(chridx)
                (corrected_fc,regions)=segmentation.do_segmentation(chrdf)
                correctedFC=correctedFC.append(corrected_fc)
                segments=segments.append(regions)
        correctedFC.reset_index(drop=True, inplace=True) # reindex data frame after concatenation
        segments.reset_index(drop=True, inplace=True)
        return (correctedFC,segments)

    def corrected_counts(nc,correctedFC,segments,grnalib,min_genes=3,ample_id='HT-29'):
        print(nc.head())
        print(correctedFC.head())
        print(segments.head())
        print(grnalib.head())
        for row in segments.itertuples():
            print (row.Index)






























































    @staticmethod
    def get_position(row):
        return round(row['startp'] + (row['endp'] - row['startp']) / 2)
    @staticmethod
    def swap_key2val(inputdict):
        return {val:key for key,val in inputdict.items()}
