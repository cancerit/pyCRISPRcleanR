
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

def runCBS(cnarr,sample_id='testSample'):
    print(R_VERSION_BUILD)
    cnarrcp=cnarr.copy()
    cnarrcp.is_copy= False
    #cnarrcp['BP'] += 1 # Convert to 1-indexed coordinates for R
    tbl = pandas2ri.py2ri(cnarrcp)
    kwargs = {'data.type':"logratio", 'sampleid':sample_id, 'presorted': True}
    # represents rx [ and rx2 [[ R brackets rx2 returns data frame
    # set seed
    base.set_seed(0xA5EED)
    cna = dnacopy.CNA(tbl.rx2('avgFC'), tbl.rx2('CHR'), tbl.rx2('BP'), **kwargs)
    smoothed_cna = dnacopy.smooth_CNA(cna)
    cnseg = dnacopy.segment(smoothed_cna, verbose=1)
    (output,data,segrows,calls)=process_output(cnseg)
    return (output,data,segrows,calls)

def process_output(cnseg):
    #ID chrom  loc.start  loc.end  num.mark  seg.mean
    output=pandas2ri.ri2py(cnseg.rx2('output'))
    #chrom  maploc  testSample
    data=pandas2ri.ri2py(cnseg.rx2('data'))
    # startRow  endRow
    segrows=pandas2ri.ri2py(cnseg.rx2('segRows'))
    calls=pandas2ri.ri2py(cnseg.rx2('call'))
    # convert coordinated back to 0 based
    #output['loc.start']-=1
    #output['loc.end']-=1
    #data['maploc']-=1
    return output,data,segrows,calls
