import sys
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.rinterface import R_VERSION_BUILD

pandas2ri.activate()

d = {'package.dependencies': 'package_dot_dependencies',
     'package_dependencies': 'package_uscore_dependencies'}

base = importr('base', robject_translations=d)
dnacopy = importr("DNAcopy", robject_translations=d)


def runCBS(cnarr, sample_id='mysample', fc_col='avgFC'):
    """
         rub CBS using DNAcopy in R
         cnseg is a final R dataframe contains 'output [#ID chrom  loc.start  loc.end
         num.mark  seg.mean], segRows[# startRow  endRow], data[#chrom  maploc  testSample]  and calls ,
         we are only interested in segRows
         to covert R data frame to python use rx [ and rx2 [[ R brackets
         rx2 returns data frame
    """
    print("CBS analysis uing R version:{}".format(R_VERSION_BUILD))
    chr_name = cnarr['CHR'].unique()[0]
    tbl = pandas2ri.py2ri(cnarr)
    kwargs = {'data.type': "logratio", 'sampleid': sample_id, 'presorted': True}
    # set seed
    base.set_seed(0xA5EED)
    cna = dnacopy.CNA(tbl.rx2(fc_col), tbl.rx2('CHR'), tbl.rx2('BP'), **kwargs)

    smoothed_cna = dnacopy.smooth_CNA(cna)
    cnseg = dnacopy.segment(smoothed_cna, verbose=1)
    segrows = pandas2ri.ri2py(cnseg.rx2('segRows'))

    return segrows, cnseg
