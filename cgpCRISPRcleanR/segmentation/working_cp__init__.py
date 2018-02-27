"""Segmentation of copy number values."""
from __future__ import absolute_import, division, print_function
from builtins import map

import locale
import logging
import math
import os.path
import numpy as np
import pandas as pd
import tempfile

from .. import core
from . import cbs

def do_segmentation(cnarr, save_dataframe=False, rlibpath=None):
    """Infer copy number segments from the given coverage table."""
    if not len(cnarr):
        return cnarr
    filtered_cn = cnarr.copy() # apply any filters before processing , should be outside this function
    # Filter out bins with no or near-zero sequencing coverage
    # Filter by distance from rolling quantiles
    # Filter by bin weights
    if not len(filtered_cn):
        return filtered_cn
    seg_out = ""
    # Run R scripts to calculate copy number segments
    rscript = {'cbs': cbs.CBS_RSCRIPT}['cbs']
    filtered_cn['BP'] += 1 # Convert to 1-indexed coordinates for R
    with tempfile.NamedTemporaryFile(suffix='.txt', mode="w+t") as tmp:
        filtered_cn.to_csv(tmp, index=False, sep='\t', float_format='%.6g', mode="w+t")
        tmp.flush()
        script_strings = {
            'chr_fc': tmp.name,
            'sample_id': 'test_sample',
            'rlibpath': ('.libPaths(c("%s"))' % rlibpath if rlibpath else ''),
        }
        with core.temp_write_text(rscript % script_strings,
                                  mode='w+t') as script_fname:
            seg_out = core.call_quiet('Rscript', '--vanilla', script_fname)
    #segarr.meta = cnarr.meta.copy()
    if save_dataframe: # for parallel execution per chromosome segment
        return seg_out
    else:
        return seg_out


#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

#from rpy2.robjects.packages import importr

#base = importr('base')
# call an R function on a Pandas DataFrame
#base.summary(my_pandas_dataframe)
