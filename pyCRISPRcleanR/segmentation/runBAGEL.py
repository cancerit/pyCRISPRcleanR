# ---------------------------------
# VERSION = 0.91
# BAGEL:  Bayesian Analysis of Gene EssentaLity
# (c) Traver Hart, 02/2015.
# modified 9/2015
# modified 25/09/2018 : Shriram  Bhosle (sb43@sanger.ac.uk)
# Free to modify and redistribute with attribtuion
# python ~/scripts/bagel-for-knockout-screens-code/BAGEL.py
# -i final_test_20_09_2018_v2/02_normalised_fold_changes.tsv
# -o bagel_out.out -e ../ref_genes/essential.txt -n ../ref_genes/non_essential.txt -c 1,2,3
# ---------------------------------

import numpy as np
import scipy.stats as stats
import sys
import logging

log = logging.getLogger(__name__)


def round_to_hundredth(x):
    return np.round(x * 100) / 100.0


def bootstrap_resample(X, n=None, ):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
    Results
    -------
    returns X_resamples

    adapted from
    Dated 7 Oct 2013
    http://nbviewer.ipython.org/gist/aflaxman/6871948
    """

    if n is None:
        n = len(X)

    # added by sb43 to achieve randomization in parallel steps
    np.random.seed(None)
    resample_i = np.floor(np.random.rand(n) * len(X)).astype(int)
    X_resample = X[resample_i]
    return X_resample


def run(bf, fc, gene_idx, genes_array, coreEss, nonEss, start, stop, thread):

    """
    :param bf: bayes factor dictionary
    :param fc: foldchanges
    :param gene_idx: gene index
    :param genes_array: gene array
    :param coreEss: essential genes
    :param nonEss: non essential genes
    :param start: bootstrap iteration window start
    :param stop: bootstrap iteration window stop
    :param thread: thread number
    :return: bayes factor dictionary
    """
    log.info("Iter", "TrainEss", "TrainNon", "TestSet")
    FC_THRESH = 2 ** -7
    sys.stdout.flush()
    count = 0
    for loop in (range(start, stop)):
        ##
        # BOOTSTRAP ITERATIONS
        #
        # bootstrap resample from gene list to get the training set
        #
        gene_train_idx = bootstrap_resample(gene_idx)
        #
        # test set for this iteration is everything not selected in bootstrap resampled training set
        #
        gene_test_idx = np.setxor1d(gene_idx, gene_train_idx)
        #
        # define essential and nonessential training sets:  arrays of indexes
        #
        count += 1
        train_ess = np.where(np.in1d(genes_array[gene_train_idx], coreEss))[0]
        train_non = np.where(np.in1d(genes_array[gene_train_idx], nonEss))[0]
        # added condition to report output after num of iterations
        if loop % 10 == 0:
            log.info("Thread:{} ITER:{},Training set Ess:{}, Non_ess:{}, test set:{}".format(thread, str(loop),
                                                                                             len(train_ess),
                                                                                             len(train_non),
                                                                                             len(gene_test_idx)))
        #
        # define ess_train: vector of observed fold changes of essential genes in training set
        #
        ess_train_fc_list_of_lists = [fc[x] for x in genes_array[gene_train_idx[train_ess]]]
        ess_train_fc_flat_list = [obs for sublist in ess_train_fc_list_of_lists for obs in sublist]
        #
        # define non_train vector of observed fold changes of nonessential genes in training set
        #
        non_train_fc_list_of_lists = [fc[x] for x in genes_array[gene_train_idx[train_non]]]
        non_train_fc_flat_list = [obs for sublist in non_train_fc_list_of_lists for obs in sublist]
        #
        # calculate empirical fold change distributions for both
        #
        kess = stats.gaussian_kde(ess_train_fc_flat_list)
        knon = stats.gaussian_kde(non_train_fc_flat_list)
        #
        # define empirical upper and lower bounds within which to calculate BF = f(fold change)
        #
        x = np.arange(-10, 2, 0.01)
        nonfitx = knon.evaluate(x)
        # define lower bound empirical fold change threshold:  minimum FC where knon is above threshold
        f = np.where(nonfitx > FC_THRESH)
        xmin = round_to_hundredth(min(x[f]))
        # define upper bound empirical fold change threshold:  minimum value of log2(ess/non)
        subx = np.arange(xmin, max(x[f]), 0.01)
        logratio_sample = np.log2(kess.evaluate(subx) / knon.evaluate(subx))
        f = np.where(logratio_sample == logratio_sample.min())
        xmax = round_to_hundredth(subx[f])
        #
        # round foldchanges to nearest 0.01
        # precalculate logratios and build lookup table (for speed)
        #
        logratio_lookup = {}
        for i in np.arange(xmin, xmax + 0.01, 0.01):
            logratio_lookup[round(i * 100)] = np.log2(kess.evaluate(i) / knon.evaluate(i))
        #
        # calculate BFs from lookup table for withheld test set
        #
        for g in genes_array[gene_test_idx]:
            foldchanges = np.array(fc[g])
            foldchanges[foldchanges < xmin] = xmin
            foldchanges[foldchanges > xmax] = xmax
            # bayes_factor = sum([logratio_lookup[round(x * 100)] for x in foldchanges])
            # sb43 modified as per fi1's recommendations
            bayes_factor = np.mean([logratio_lookup[round(x * 100)] for x in foldchanges])
            bf[g].append(bayes_factor)
    return bf
