from multiprocessing import Pool
from . import cbs
from . import runBAGEL
import numpy as np
import logging

log = logging.getLogger(__name__)

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


def run_parallel_bagel(foldchangefile, column_list, Ess, nonEss, cpus, NUM_BOOTSTRAPS=1000):
    # DEFINE REFERENCE SETS
    coreEss = np.array(Ess)
    log.info("Number of reference essentials: {} ".format(len(coreEss)))
    nonEss = np.array(nonEss)
    log.info("Number of reference nonessentials: {} ".format(len(nonEss)))
    bf, bf_dict, fc, gene_idx, genes_array = _prepare_data(foldchangefile, column_list)
    boot_range = list(map(int, np.linspace(0, NUM_BOOTSTRAPS, cpus + 1).tolist()))
    with Pool(cpus) as (pool):
        bf_pool = list(pool.map(_process_df,
                                ((bf, fc, gene_idx, genes_array, coreEss, nonEss,
                                  boot_range, boot_iter) for boot_iter in range(0, cpus, 1))))
    for bf in bf_pool:
        for gene, baf in bf.items():
            if baf:
                bf_dict[gene] += baf
    log.info("Completed bagel analysis")
    return bf_dict


def _process_df(args):
    """ do things in parallel"""
    return _run_bagel(*args)


def _run_bagel(bf, fc, gene_idx, genes_array, coreEss, nonEss, boot_range, boot_iter):
    start = boot_range[boot_iter]
    stop = boot_range[boot_iter + 1]
    log.info("Running bagel iteration between: {} and {} ".format(start, stop))
    bf_res = runBAGEL.run(bf, fc, gene_idx, genes_array, coreEss, nonEss, start, stop, boot_iter)
    return bf_res


def _prepare_data(foldchangefile, column_list):
    # LOAD FOLDCHANGES
    genes = {}
    fc = {}
    fin = open(foldchangefile)
    skipfields = fin.readline().rstrip().split('\t')
    for i in column_list:
        print("Using column:{}".format(skipfields[i + 1]))
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
    log.info("Number of unique genes: {} ".format(len(genes)))

    # INITIALIZE BFS
    #
    bf = {}
    bf_dict = {}
    for g in genes_array:
        bf[g] = []
        bf_dict[g] = []
    return bf, bf_dict, fc, gene_idx, genes_array
