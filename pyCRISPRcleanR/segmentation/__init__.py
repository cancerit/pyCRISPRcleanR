from multiprocessing import Pool
from . import cbs
import numpy as np

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


def run_bagel(foldchangefile, column_list, Ess, nonEss, cpus, BOOTSTRAPS=1000):
    bf_dict = {}
    bf, fc, gene_idx, genes_array, coreEss, nonEss = _prepare_data(foldchangefile, column_list, noEss, nonEss)
    with Pool(cpus) as (pool):
        result = list(pool.map(_ds, ((bf, fc, gene_idx, genes_array, coreEss, nonEss, NUM_BOOTSTRAPS=boot_iter) for
                                     boot_iter in range(BOOTSTRAPS))))
        for result_dict in result:
            bf_dict.update(result_dict)
    return bf_dict


def _prepare_data(foldchangefile, column_list, coreEss, nonEss):
    # LOAD FOLDCHANGES
    genes = {}
    fc = {}
    fin = open(foldchangefile)
    skipfields = fin.readline().rstrip().split('\t')
    for i in column_list:
        print("Using column:" + skipfields[i + 1])
    for line in fin:
        fields = line.rstrip().split('\t')
        gsym = fields[1]
        genes[gsym] = 1
        if gsym not in fc:
            fc[gsym] = []  # initialize dict entry as a list
        for i in column_list:
            # per user docs, GENE is column 0, first data column is col 1.
            fc[gsym].append(float(fields[i + 1]))
    fin.close()  # sb43 added filehandle closure
    genes_array = np.array(list(genes.keys()))
    gene_idx = np.arange(len(genes))
    # print "Number of gRNA loaded:  " + str( len(genes_array) )
    print("Number of unique genes:  " + str(len(genes)))

    # DEFINE REFERENCE SETS
    coreEss = np.array(coreEss)
    print("Number of reference essentials: " + str(len(coreEss)))
    nonEss = np.array(nonEss)
    print("Number of reference nonessentials: " + str(len(nonEss)))
    #
    # INITIALIZE BFS
    #
    bf = {}
    for g in genes_array:
        bf[g] = []
    return bf, fc, gene_idx, genes_array, coreEss, nonEss
