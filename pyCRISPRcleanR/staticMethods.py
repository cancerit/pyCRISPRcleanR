import sys
import os
import tarfile
import pandas as pd
import numpy as np
import logging
from pyCRISPRcleanR.plots import PlotData as PLT
from . import segmentation

log = logging.getLogger(__name__)


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

    # ------------------------------Analysis methods---------------------------------
    @staticmethod
    def combine_count_n_library(countfile, libfile, plot_flag=None, outdir='./'):
        """
            Combine counts and library file data based on union of indexes
        """
        try:
            counts = pd.read_csv(countfile, compression='infer', sep="\t", index_col='sgRNA')
            libdata = pd.read_csv(libfile, compression='infer', sep="\t", index_col='sgRNA')
        except ValueError:
            sys.exit('Invalid index column name, please check input format for count and library file')

        if {'gene'}.issubset(counts.columns):
            counts.drop(['gene'], axis=1, inplace=True, errors='raise')
        else:
            sys.exit('counts data does not contain required column:[gene]')
        # plot raw counts
        if plot_flag:
            PLT.box_plot_r(counts.iloc[:, 1:], title="Raw sgRNA counts", saveto=outdir + '/raw_counts',
                           ylabel='Raw Counts', xlabel='Sample Names')
            PLT.box_plot_ly(counts.iloc[:, 1:], title="Raw sgRNA counts", saveto=outdir + '/raw_counts_plotly',
                            ylabel='Raw Counts', xlabel='Sample Names')
            log.info("Plotted raw counts.....")

        if {'gene', 'chr', 'start', 'end'}.issubset(libdata.columns):
            libdata = libdata.get(['gene', 'chr', 'start', 'end'])
        else:
            sys.exit('Library data does not contain required columns:[gene, chr, start, end]')

        cldf = pd.concat([counts, libdata], axis=1, join='inner')
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
    def get_norm_count_n_fold_changes(cldf, controls, plot_flag=None, outdir='./'):
        """
            Calculate normalised count and avarage fold change using raw count data for
            control and sample, add  it to main data frame
            Drop control fold change column so that average fold change can be calculated on
            sample fold changes only
            Sort dataframe using chr and start position of gRNA
        """
        normed = cldf.iloc[:, 0:cldf.columns.get_loc('gene')].div(
            cldf.iloc[:, 0:cldf.columns.get_loc('gene')].agg('sum')) * 10e6
        # plot normalised counts
        if normed.empty:
            sys.exit('Normalized data frame is empty check if required columns are present')

        elif plot_flag:
            PLT.box_plot_r(normed, title="Normalised sgRNA counts", saveto=outdir + '/normalised_counts',
                           ylabel='Normalised Counts',
                           xlabel='Sample Names')
            PLT.box_plot_ly(normed, title="Normalised sgRNA counts", saveto=outdir + '/normalised_counts',
                            ylabel='Normalised Counts',
                            xlabel='Sample Names')

        fc = normed.apply(lambda x: np.log2((x + 0.5) / (normed.iloc[:, 0:controls].mean(axis=1) + 0.5)))
        # drop control columns
        fc.drop(fc.columns[0:controls], axis=1, inplace=True)

        if fc.empty:
            sys.exit('Foldchange data frame is empty')
            # plot fold Changes
        elif plot_flag:
            PLT.box_plot_r(fc, title="Fold Changes sgRNA", saveto=outdir + '/fold_changes',
                           ylabel='Fold Changes',
                           xlabel='Sample Names')
            PLT.box_plot_ly(fc, title="Fold Changes sgRNA", saveto=outdir + '/fold_changes_plotly',
                            ylabel='Fold Changes',
                            xlabel='Sample Names')

        no_rep = len(fc.columns)

        cldf = cldf.join(normed, rsuffix='_nc')
        normed.insert(0, 'gene', cldf['gene'])
        StaticMthods._print_df(normed, outdir + "/normalised_counts.tsv")

        cldf['avgFC'] = fc.mean(axis=1)
        fc['avgFC'] = cldf['avgFC']
        fc.insert(0, 'gene', cldf['gene'])
        StaticMthods._print_df(fc, outdir + "/normalised_fold_changes.tsv")

        cldf['BP'] = round(cldf['start'] + (cldf['end'] - cldf['start']) / 2).astype(int)
        cldf.sort_values(by=['chr', 'start'], ascending=True, inplace=True)
        return cldf, no_rep

    @staticmethod
    def run_cbs(cldf, cpus, sample, fc_col='avgFC'):
        """
            Runs CBS algorithm from DNAcopy and returns a dictionay
            of per chr raw dataframe and cbs segments

        """
        cbs_dict = segmentation.do_segmentation(cldf, cpus, sample, fc_col=fc_col)
        return cbs_dict

    @staticmethod
    def process_segments(cbs_dict, ignored_genes, min_genes, controls, no_rep, outdir='./'):
        """
           Process CBS derived copy number segment to get Correted foldchanges and
           reverted count data
        """
        corrected_count_list = []
        chrdata_list = []
        for chr, (cnarr, segrows, cnseg) in cbs_dict.items():
            print("Correcting counts on Chr :{} : sgRNA:{} segments{}".format(chr, cnarr.shape, segrows.shape))
            cnarr.is_copy = False
            cnarr['correction'] = 0
            cnarr['correctedFC'] = cnarr.avgFC
            n_genes_in_seg = 0
            reverted_counts = cnarr.iloc[:, cnarr.columns.get_loc('end') + controls + 1:
                                        cnarr.columns.get_loc('avgFC')]

            for segment in segrows.itertuples():
                idxs = list(range(segment.startRow - 1, segment.endRow))
                included_genes = cnarr.gene.iloc[idxs].unique()
                n_genes_in_seg = len(set(included_genes) - set(ignored_genes))
                if n_genes_in_seg >= min_genes:
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
        # get control counts to and join with corrected_counts for printing
        nc_control_count = alldata.iloc[:, alldata.columns.get_loc('end') + 1:
                                    alldata.columns.get_loc('end') + controls + 1]
        corrected_count = nc_control_count.join(corrected_count)
        corrected_count = corrected_count.rename(columns=lambda x: str(x)[:-3])
        # calculate corrected fold changes
        corrected_fc = corrected_count.apply(
            lambda x: np.log2((x + 0.5) / (corrected_count.iloc[:, 0:controls].mean(axis=1) + 0.5)))
        corrected_fc.drop(corrected_fc.columns[0:controls], axis=1, inplace=True)
        corrected_fc['avgFC'] = corrected_fc.mean(axis=1)
        alldata = alldata.join(corrected_count, rsuffix='_cc')
        alldata = alldata.join(corrected_fc, rsuffix='_cf')

        # add gene names before writing to a file
        corrected_fc.insert(0, 'gene', alldata['gene'])
        StaticMthods._print_df(corrected_fc, outdir + "/crispr_cleanr_fold_changes.tsv")
        # add gene names before writing to a file
        corrected_count.insert(0, 'gene', alldata['gene'])
        StaticMthods._print_df(corrected_count, outdir + "/crispr_cleanr_corrected_counts.tsv")
        return alldata

    @staticmethod
    def _correct_counts(segdata, controls, no_rep):
        """
          correct count based on segmentation output
        """
        reverted = pd.DataFrame()
        # average fold change
        # mean normalised control count
        nc = segdata.iloc[:, segdata.columns.get_loc('end') + 1:segdata.columns.get_loc('end') + controls + 1]
        c = nc.mean(axis=1)
        n = segdata.correctedFC
        reverted['revc'] = c * (pow(2, n))
        normed_num = segdata.iloc[:, segdata.columns.get_loc('end') + controls + 1:
                                segdata.columns.get_loc('avgFC')]
        normed_num += 1
        proportions = normed_num.div(normed_num.agg('sum', axis=1), axis=0)
        reverted = reverted * no_rep
        reverted = proportions.mul(reverted.revc, axis=0)
        return reverted

    @staticmethod
    def _print_df(mydf, out_file):
        log.info("Writing out file  .....:{}".format(out_file))
        mydf.to_csv(out_file, sep='\t', mode='w', header=True,
                    index=True, index_label='sgRNA', doublequote=False)
