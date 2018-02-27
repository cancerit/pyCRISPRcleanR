import sys
import os
import tarfile
import tempfile
import pandas as pd
from subprocess import Popen, PIPE, STDOUT
import shlex
import re
import time
import logging.config
#os.environ["LD_LIBRARY_PATH"] = "/software/R-3.3.0/lib/R/lib"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/R-3.3.0/lib/R/li

from . import segmentation

from beautifultable import BeautifulTable

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
        for chridx in fc.CHR.unique():
            if chridx == "1":
                chrdf=fc[fc.CHR == chridx ]
                print(chridx)
                segdata=segmentation.do_segmentation(chrdf)
                #print(segdata.shape)
        return None


#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

#from rpy2.robjects.packages import importr

#base = importr('base')
# call an R function on a Pandas DataFrame
#base.summary(my_pandas_dataframe)





    @staticmethod
    def get_position(row):
        return round(row['startp'] + (row['endp'] - row['startp']) / 2)
    @staticmethod
    def swap_key2val(inputdict):
        return {val:key for key,val in inputdict.items()}
