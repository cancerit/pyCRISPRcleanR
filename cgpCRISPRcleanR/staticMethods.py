import sys
import os
import tarfile
import pandas as pd
import numpy as np
import logging.config
# os.environ["LD_LIBRARY_PATH"] = "/software/R-3.3.0/lib/R/lib"
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/R-3.3.0/lib/R/li

from . import segmentation

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
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

    @staticmethod
    def create_dict(inputlist):
        chrdict = {}
        chr_count = 1
        for chr in inputlist:
            chrdict[chr] = chr_count
            chr_count += 1
        return chrdict

    # ------------------------------Analysis methods---------------------------------
    @staticmethod
    def combine_count_n_library(countfile, libfile):
        """
            Combine counts and library file data based on union of indexes
        """
        counts = pd.read_csv(countfile, compression='infer', sep="\t", index_col='sgRNA')
        libdata = pd.read_csv(libfile, compression='infer', sep="\t", index_col='sgRNA')
        cldf = pd.concat([counts, libdata], axis=1, join='inner')
        return cldf

    @staticmethod
    def cleanup_data(cldf):
        """
            perform cleanup: change row names and drop duplicate columns
            This method may not be required once input format for library and
            counts table is standardized
        """
        col_to_drop = ['gene', 'EXONE', 'CODE', 'STRAND']
        rename_col_dict = {'CHRM': 'CHR', 'STARTpos': 'startp', 'ENDpos': 'endp', 'GENES': 'gene'}
        cldf.drop(col_to_drop, axis=1, inplace=True, errors='raise')
        cldf = cldf.rename(columns=rename_col_dict)
        return cldf

    @staticmethod
    def filter_data(cldf, controls, min_read_count):
        """
            filter data frame for minimum read counts cutoff in control sample
            Add any other filter criteria in future....
        """
        cldf.drop(cldf[cldf.iloc[:, 0:controls].mean(axis=1) < min_read_count].index, inplace=True)
        return cldf

    @staticmethod
    def get_norm_count_n_fold_changes(cldf, controls):
        """
            Calculate normalised count and avarage fold change using raw count data for
            control and sample, add  it to main data frame
            Drop control fold change column so that average fold change can be calculated on
            sample fold changes only
            Sort dataframe using chr and start position of gRNA
        """
        normed = cldf.iloc[:, 0:cldf.columns.get_loc('gene')].div(
            cldf.iloc[:, 0:cldf.columns.get_loc('gene')].agg('sum')) * 10e6
        if normed.empty:
            print(normed.head())
            sys.exit('Normalized data frame is empty check if required columns are present')
        fc = normed.apply(lambda x: np.log2((x + 0.5) / (normed.iloc[:, 0:controls].mean(axis=1) + 0.5)))
        fc.drop(fc.columns[0:controls], axis=1, inplace=True)
        if fc.empty:
            print(fc.head())
            sys.exit('Foldchange data frame is empty')
        no_rep = len(fc.columns)
        cldf = cldf.join(normed, rsuffix='_nc')
        cldf['avgFC'] = fc.mean(axis=1)
        cldf['BP'] = round(cldf['startp'] + (cldf['endp'] - cldf['startp']) / 2).astype(int)
        cldf.sort_values(by=['CHR', 'startp'], ascending=True, inplace=True)
        return cldf, no_rep

    @staticmethod
    def run_cbs(cldf, cpus, sample):
        """
            Runs CBS algorithm from DNAcopy and returns a dictionay
            of per chr raw dataframe and cbs segments

        """
        cbs_dict = segmentation.do_segmentation(cldf, cpus, sample)
        return cbs_dict

    @staticmethod
    def process_segments(cbs_dict, ignored_genes, min_genes, controls, no_rep):
        """
           Process CBS derived copy number segment to get Correted foldchanges and
           reverted count data
        """
        corrected_count_list = []
        chrdata_list = []
        for chr, (cnarr, segrows) in cbs_dict.items():
            print("Correcting counts on Chr :{} : sgRNA:{} segments{}".format(chr, cnarr.shape, segrows.shape))
            cnarr.is_copy = False
            cnarr['correction'] = 0
            cnarr['correctedFC'] = cnarr.avgFC
            reverted_counts = cnarr.iloc[:, cnarr.columns.get_loc('endp') +
                            controls + 1:cnarr.columns.get_loc('avgFC')]
            n_gene_in_seg = 0
            for segment in segrows.itertuples():
                idxs = list(range(segment.startRow - 1, segment.endRow))
                included_genes = cnarr.gene.iloc[idxs].unique()
                n_gene_in_seg = len(set(included_genes) - set(ignored_genes))
                if n_gene_in_seg >= min_genes:
                    cnarr.iloc[(idxs, cnarr.columns.get_loc('correctedFC'))] = \
                        cnarr.avgFC.iloc[idxs] - cnarr.avgFC.iloc[idxs].mean()
                    cnarr.iloc[(idxs, cnarr.columns.get_loc('correction'))] = -np.sign(
                        cnarr.correctedFC.iloc[idxs].mean())
                    reverted = StaticMthods._correct_counts(cnarr.iloc[idxs], controls, no_rep)
                    reverted_counts.iloc[idxs] = reverted
            corrected_count_list.append(reverted_counts)
            chrdata_list.append(cnarr)
        corrected_count = pd.concat(corrected_count_list)
        alldata = pd.concat(chrdata_list)
        alldata = alldata.join(corrected_count, rsuffix='_rev')
        return alldata

    @staticmethod
    def _correct_counts(segdata, controls, no_rep):
        """
          correct count based on segmentation output
        """
        reverted = pd.DataFrame()
        # average fold change
        # mean normalised control count
        nc = segdata.iloc[:, segdata.columns.get_loc('endp') + 1:segdata.columns.get_loc('endp') + controls + 1]
        c = nc.mean(axis=1)
        n = segdata.correctedFC
        reverted['revc'] = c * (pow(2, n))
        normed_num = segdata.iloc[:, segdata.columns.get_loc('endp') +
                                controls + 1:segdata.columns.get_loc('avgFC')]
        normed_num += 1
        proportions = normed_num.div(normed_num.agg('sum', axis=1), axis=0)
        reverted = reverted * no_rep
        reverted = proportions.mul(reverted.revc, axis=0)
        return reverted

# utility methods
    @staticmethod
    def get_position(row):
        return round(row['startp'] + (row['endp'] - row['startp']) / 2)

    @staticmethod
    def swap_key2val(inputdict):
        return {val: key for key, val in inputdict.items()}
