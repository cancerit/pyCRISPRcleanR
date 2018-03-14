
import numpy as np
import pandas as pd

from rpy2.robjects import r, pandas2ri
from rpy2.rinterface import R_VERSION_BUILD
pandas2ri.activate()
from rpy2.robjects.packages import importr

d = {'package.dependencies': 'package_dot_dependencies',
     'package_dependencies': 'package_uscore_dependencies'}

base = importr('base',robject_translations=d)
#base.set.seed(0xA5EED)
dnacopy=importr("DNAcopy",robject_translations=d)

def runCBS(cnarr,sample_id='mysample'):
    """
         rub CBS using DNAcopy in R
         cnseg is a final R dataframe contains 'output [#ID chrom  loc.start  loc.end  
         num.mark  seg.mean] , 
        , segRows[# startRow  endRow], data[#chrom  maploc  testSample]  and calls ,
         we are only interested in segRows
         to covert R data frame to python use rx [ and rx2 [[ R brackets 
         rx2 returns data frame
    """
    print("CBS analysis uing R version:{}".format(R_VERSION_BUILD))
    tbl = pandas2ri.py2ri(cnarr)
    kwargs = {'data.type':"logratio", 'sampleid':sample_id, 'presorted': True}
    # set seed
    base.set_seed(0xA5EED)
    cna = dnacopy.CNA(tbl.rx2('avgFC'), tbl.rx2('CHR'), tbl.rx2('BP'), **kwargs)
    smoothed_cna = dnacopy.smooth_CNA(cna)
    cnseg = dnacopy.segment(smoothed_cna, verbose=1)
    segrows=pandas2ri.ri2py(cnseg.rx2('segRows'))
    return segrows

