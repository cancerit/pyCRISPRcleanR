import sys
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import logging
import numpy as np

log = logging.getLogger(__name__)

pandas2ri.activate()
d = {'package.dependencies': 'package_dot_dependencies',
     'package_dependencies': 'package_uscore_dependencies'}
base = importr('base', robject_translations=d)
log.info(base._libPaths())
dnacopy = importr("DNAcopy", robject_translations=d)


def runCBS(cnarr, fc_col='avgFC'):
    """
         rub CBS using DNAcopy in R
         cnseg is a final R dataframe contains 'output [#ID chrom  loc.start  loc.end
         num.mark  seg.mean], segRows[# startRow  endRow], data[#chrom  maploc  testSample]  and calls ,
         we are only interested in segRows
         to covert R data frame to python use rx [ and rx2 [[ R brackets
         rx2 returns data frame
    """
    chr_name = cnarr['chr'].unique()[0]
    tbl = pandas2ri.py2ri(cnarr)
    kwargs = {'data.type': "logratio", 'presorted': True}
    # set seed
    base.set_seed(0xA5EED)
    cna = dnacopy.CNA(tbl.rx2(fc_col), tbl.rx2('chr'), tbl.rx2('BP'), **kwargs)
    
    smoothed_cna = dnacopy.smooth_CNA(cna)

    cnseg = dnacopy.segment(smoothed_cna, verbose=1)
    segrows = pandas2ri.ri2py(cnseg.rx2('segRows'))
    return segrows, cnseg
