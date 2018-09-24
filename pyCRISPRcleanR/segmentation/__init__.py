from multiprocessing import Pool
from . import cbs

"""Segmentation of loagratio fold change  values."""


def do_segmentation(cnarr, cpus, fc_col='avgFC'):
    cnseg_dict = {}
    with Pool(cpus) as (pool):
        result = list(pool.map(_ds, ((chrname, ca, fc_col) for chrname, ca in cnarr.groupby('chr'))))
        for result_dict in result:
            cnseg_dict.update(result_dict)
    return cnseg_dict


def _ds(args):
    """Wrapper for parallel map"""
    return _do_segmentation(*args)


def _do_segmentation(chrname, ca, fc_col):
    """
    :rtype: result dictionary
    """
    result_dict = {}
    print(('Performing CBS on chr:{}').format(chrname))
    segrows, cnseg = cbs.runCBS(ca, fc_col=fc_col)
    result_dict[chrname] = [ca, segrows, cnseg]
    return result_dict
