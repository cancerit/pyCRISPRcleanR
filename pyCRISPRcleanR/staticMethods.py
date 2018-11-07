import sys
import os
import tarfile
import json
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import numpy as np
import logging
from pyCRISPRcleanR.plots import PlotData as PLT
from . import segmentation

log = logging.getLogger(__name__)

MAGECK_CMD = "mageck test --count-table {} --control-id {} --treatment-id {} --output-prefix {} --norm-method {}"
SIGNATURE_FILES = ("essential", "non_essential", "dna_replication", "rna_polymerase",
                   "proteasome", "ribosomal_proteins", "spliceosome")
CONTROL_SAMPLES = 'NA'
TREATMENT_SAMPLES = 'NA'

RESULTS_FILE = 'results'

RAW_BOXPLOT = "/01_raw_counts_boxplot"
RAW_HIST = "/01_raw_counts_histogram"
RAW_MATRIX = "/01_raw_counts_correlation_matrix"
NORM_BOXPLOT = "/02_normalised_counts_boxplot"
NORM_HIST = "/02_normalised_counts_histogram"
NORM_MATRIX = "/02_normalised_counts_correlation_matrix"
FC_BOXPLOT = "/03_fold_changes_boxplot"
FC_HIST = "/03_fold_changes_histogram"
FC_MATRIX = "/03_fold_changes_correlation_matrix"
CC_BOXPLOT = "/07_CRISPRcleanR_corrected_count_boxplot"
CC_HIST = "/07_CRISPRcleanR_corrected_count_histogram"
CC_MATRIX = "/07_CRISPRcleanR_corrected_count_correlation_matrix"
CFC_BOXPLOT = "/08_CRISPRcleanR_corrected_fold_changes_boxplot"
CFC_HIST = "/08_CRISPRcleanR_corrected_fold_changes_histogram"
CFC_MATRIX = "/08_CRISPRcleanR_corrected_fold_changes_correlation_matrix"

NC_TSV = "/01_normalised_counts.tsv"
NFC_TSV = "/02_normalised_fold_changes.tsv"
CFC_TSV = "/04_crispr_cleanr_fold_changes.tsv"
CC_TSV = "/03_crispr_cleanr_corrected_counts.tsv"


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
    def combine_count_n_library(countfile, libfile, outdir='./'):
        """
         Combine counts and library file data based on union of indexes
        :param countfile:
        :param libfile:
        :param outdir:
        :return:
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
        PLT.box_plot_ly(counts, title="Raw sgRNA counts", saveto=outdir + RAW_BOXPLOT,
                        ylabel='Raw Counts', xlabel='Sample Names')
        PLT.histogram_ly(counts, title="Raw sgRNA counts", saveto=outdir + RAW_HIST,
                         ylabel='Raw Counts',
                         xlabel='sgRNAbins')
        PLT.correlation_plot_ly(counts, title="Correlation matrix: raw sgRNA counts",
                                saveto=outdir + RAW_MATRIX,
                                ylabel='Raw Counts',
                                xlabel='Sample Names')
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
    def get_norm_count_n_fold_changes(cldf, controls, outdir='./'):
        """
            Calculate normalised count and avarage fold change using raw count data for
            control and sample, add  it to main data frame
            Drop control fold change column so that average fold change can be calculated on
            sample fold changes only
            Sort dataframe using chr and start position of gRNA
        :param cldf: counts and annotation library df
        :param controls:
        :param outdir:
        :return:
        """
        global CONTROL_SAMPLES
        global TREATMENT_SAMPLES
        global NORM_HIST, NORM_BOXPLOT, NORM_MATRIX, FC_MATRIX, FC_HIST, FC_BOXPLOT, NC_TSV, NFC_TSV

        normed = cldf.iloc[:, 0:cldf.columns.get_loc('gene')].div(
            cldf.iloc[:, 0:cldf.columns.get_loc('gene')].agg('sum')) * 10e6
        # plot normalised counts
        if normed.empty:
            sys.exit('Normalized data frame is empty check if required columns are present')

        PLT.box_plot_ly(normed, title="Normalised sgRNA counts",
                        saveto=outdir + NORM_BOXPLOT,
                        ylabel='Normalised Counts',
                        xlabel='Sample Names')

        PLT.histogram_ly(normed, title="Normalised sgRNA counts",
                         saveto=outdir + NORM_HIST,
                         ylabel='Normalised Counts',
                         xlabel='sgRNAbins')
        PLT.correlation_plot_ly(normed, title="Correlation matrix: normalised sgRNA counts",
                                saveto=outdir + NORM_MATRIX,
                                ylabel='Normalised Counts',
                                xlabel='Sample Names')

        fc = normed.apply(lambda x: np.log2((x + 0.5) / (normed.iloc[:, 0:controls].mean(axis=1) + 0.5)))
        # drop control columns
        fc.drop(fc.columns[0:controls], axis=1, inplace=True)

        if fc.empty:
            sys.exit('Foldchange data frame is empty')

            # plot fold Changes

        PLT.box_plot_ly(fc, title="Fold Changes sgRNA", saveto=outdir + FC_BOXPLOT,
                        ylabel='Fold Changes',
                        xlabel='Sample Names')
        PLT.histogram_ly(fc, title="Fold changes sgRNA", saveto=outdir + FC_HIST,
                         ylabel='Fold changes',
                         xlabel='sgRNAbins')
        PLT.correlation_plot_ly(fc, title="Correlation matrix: Fold changes sgRNA",
                                saveto=outdir + FC_MATRIX,
                                ylabel='Fold changes',
                                xlabel='Sample Names')
        no_rep = len(fc.columns)

        cldf = cldf.join(normed, rsuffix='_nc')
        normed.insert(0, 'gene', cldf['gene'])
        StaticMthods._print_df(normed, outdir + NC_TSV)

        # MAGECK  specific prms.....
        CONTROL_SAMPLES = ",".join(normed.columns[1:controls + 1].values.tolist())
        TREATMENT_SAMPLES = ",".join(normed.columns[controls + 1:].values.tolist())
        # end MAGECK specific prms

        cldf['avgFC'] = fc.mean(axis=1)
        fc['avgFC'] = cldf['avgFC']
        fc.insert(0, 'gene', cldf['gene'])
        StaticMthods._print_df(fc, outdir + NFC_TSV)
        # mean fold changes grouped by gene
        geneFC = fc.groupby(['gene'])['avgFC'].mean()
        sgRNAFC = fc['avgFC']
        cldf['BP'] = round(cldf['start'] + (cldf['end'] - cldf['start']) / 2).astype(int)
        cldf.sort_values(by=['chr', 'start'], ascending=True, inplace=True)

        return cldf, no_rep, outdir + NC_TSV, geneFC, sgRNAFC, fc, outdir + NFC_TSV

    @staticmethod
    def run_cbs(cldf, cpus, fc_col='avgFC'):
        """
            Runs CBS algorithm from DNAcopy and returns a dictionay
            of per chr raw dataframe and cbs segments

        """
        cbs_dict = segmentation.do_segmentation(cldf, cpus, fc_col=fc_col)
        return cbs_dict

    @staticmethod
    def run_bagel(fc_file, Ess, nonEss, cpus, column_list=[], NUM_BOOTSTRAPS=1000, outfilename='./bagel_out'):
        """
        :param fc_file: file containing crispr foldchange data as defined in BAGEL documentation
        :param Ess: essential genes
        :param nonEss: non essential genes
        :param cpus: num of cpus
        :param column_list: columns to analyse [see BAGEL documentation]
        :param NUM_BOOTSTRAPS: number of bootstrap iterations
        :param outfilename: output file
        :return:  return on success
        """
        bf = segmentation.run_parallel_bagel(fc_file, column_list, Ess, nonEss,
                                             cpus, NUM_BOOTSTRAPS=NUM_BOOTSTRAPS)
        fout = open(outfilename, 'w')
        fout.write('GENE\tBF\tSTD\tNumObs\n')
        for g in sorted(bf.keys()):
            if bf[g]:  # sb43 added to avoid empty array error
                num_obs = len(bf[g])
                bf_mean = np.mean(bf[g])
                bf_std = np.std(bf[g])
                fout.write('{0:s}\t{1:4.3f}\t{2:4.3f}\t{3:d}\n'.format(g, bf_mean, bf_std, num_obs))
        fout.close()
        return 0

    @staticmethod
    def process_segments(cbs_dict, ignored_genes, min_genes, controls, no_rep, outdir='./'):
        """
        :Process CBS derived copy number segment to get Correted foldchanges and reverse transformed counts
        :param cbs_dict:
        :param ignored_genes:
        :param min_genes:
        :param controls:
        :param no_rep:
        :param outdir:
        :return:
        """
        corrected_count_list = []
        chrdata_list = []
        global CC_TSV, CC_MATRIX, CC_HIST, CC_BOXPLOT, CFC_TSV, CFC_MATRIX, CFC_HIST, CFC_BOXPLOT
        for chr, (cnarr, segrows, cnseg) in cbs_dict.items():
            log.info("Correcting counts on Chr :{} : sgRNA:{} segments{}".format(chr, cnarr.shape, segrows.shape))
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

        PLT.box_plot_ly(corrected_count, title="CRISPRcleanR corrected count sgRNA",
                        saveto=outdir + CC_BOXPLOT,
                        ylabel='count',
                        xlabel='Sample Names')
        PLT.histogram_ly(corrected_count, title="CRISPRcleanR corrected count sgRNA",
                         saveto=outdir + CC_HIST,
                         ylabel='count', xlabel='sgRNAbins')
        PLT.correlation_plot_ly(corrected_count,
                                title="Correlation matrix: CRISPRcleanR corrected count sgRNA",
                                saveto=outdir + CC_MATRIX,
                                ylabel='count', xlabel='Sample Names')

        PLT.box_plot_ly(corrected_fc, title="CRISPRcleanR corrected Fold Changes sgRNA",
                        saveto=outdir + CFC_BOXPLOT,
                        ylabel='Fold Changes', xlabel='Sample Names')
        PLT.histogram_ly(corrected_fc, title="CRISPRcleanR Fold changes sgRNA",
                         saveto=outdir + CFC_HIST,
                         ylabel='Fold changes', xlabel='sgRNAbins')
        PLT.correlation_plot_ly(corrected_fc,
                                title="Correlation matrix: CRISPRcleanR corrected Fold changes sgRNA",
                                saveto=outdir + CFC_MATRIX,
                                ylabel='Fold changes', xlabel='Sample Names')

        corrected_fc['avgFC'] = corrected_fc.mean(axis=1)
        alldata = alldata.join(corrected_count, rsuffix='_cc')
        alldata = alldata.join(corrected_fc, rsuffix='_cf')
        # add gene names before writing to a file
        corrected_fc.insert(0, 'gene', alldata['gene'])
        StaticMthods._print_df(corrected_fc, outdir + CFC_TSV)
        # add gene names before writing to a file
        corrected_count.insert(0, 'gene', alldata['gene'])
        StaticMthods._print_df(corrected_count, outdir + CC_TSV)

        return alldata, outdir + CC_TSV, corrected_fc, outdir + CFC_TSV

    @staticmethod
    def _correct_counts(segdata, controls, no_rep):
        """
        :correct count based on segmentation output
        :param segdata:
        :param controls:
        :param no_rep:
        :return: reverted: counts
        """
        reverted = pd.DataFrame()
        # average fold change
        # mean normalised control count
        nc = segdata.iloc[:, segdata.columns.get_loc('end') + 1:segdata.columns.get_loc('end') + controls + 1]
        c = nc.mean(axis=1)
        n = segdata.correctedFC
        reverted['revc'] = c * (pow(2, n))
        normed_num = segdata.iloc[:, segdata.columns.get_loc('end') + controls +
                                  1:segdata.columns.get_loc('avgFC')]
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

    @staticmethod
    def run_mageck(count_file, outfile_prefix="./mageckOut/mageck", ):
        """
        :param norm_count_file:
        :param corrected_count_file:
        :param outdir:
        :return:
        """
        global CONTROL_SAMPLES
        global TREATMENT_SAMPLES

        cmd = MAGECK_CMD.format(count_file, CONTROL_SAMPLES, TREATMENT_SAMPLES,
                                outfile_prefix, 'none')
        StaticMthods._run_command(cmd)
        return outfile_prefix + '.gene_summary.txt'

    @staticmethod
    def _run_command(cmd):
        """
        : runs external command
        :param cmd:
        :return: command output
        """
        """ runs command in a shell, returns stdout and exit code"""
        if not len(cmd):
            raise ValueError("Must supply at least one argument")
        try:
            # To capture standard error in the result, use stderr=subprocess.STDOUT:
            cmd_obj = Popen(cmd, stdin=None, stdout=PIPE, stderr=PIPE,
                            shell=True, universal_newlines=True, bufsize=-1,
                            close_fds=True, executable='/bin/bash')
            log.info("running command:{}".format(cmd))
            (out, error) = cmd_obj.communicate()
            exit_code = cmd_obj.returncode
            if (exit_code == 0):
                log.info("mageck run successfully")
            else:
                log.debug("Error: mageck exited with non zero exit status, please check log file more details")
                log.error("OUT:{}:Error:{}:Exit:{}".format(out, error, exit_code))
            return 0
        except OSError as oe:
            log.error("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))
            sys.exit("Unable to run command:{} Error:{}".format(cmd, oe.args[0]))

    @staticmethod
    def load_signature_files(sig_dir_path, df):
        """
        :param sig_dir_path:
        :return: signature_dict
        """
        signature_dict = {}
        signature = None
        try:
            for signature in SIGNATURE_FILES:
                with open(sig_dir_path + '/' + signature + '.txt') as f:
                    gene_list = f.read().splitlines()
                    sgRNA_list = df[df.gene.isin(gene_list)].index.tolist()
                    signature_dict[signature + '_sgRNAs'] = sgRNA_list
                    signature_dict[signature + '_genes'] = gene_list
        except IOError:
            log.error("Unable to load signature file for:" + signature)
        return signature_dict

    @staticmethod
    def get_data_for_density_plot(df, essential, non_essential):
        """
        :param df: pandas final df
        :param essential: list of essential gRNAs
        :param non_essential: list of non_essential gRNAs
        :return:
        """
        essential_df = df[df.index.isin(essential)]
        non_essential_df = df[df.index.isin(non_essential)]
        other_df = df[~df.index.isin(essential + non_essential)]
        return essential_df, non_essential_df, other_df

    @staticmethod
    def get_obs_predictions(df_data, positive, negative):
        """
        :param df_data: numpy array
        :param positive:  positive gRNA set
        :param negative: negative gRNA set
        :return:
        """
        df = df_data.to_frame()
        df = df[df.index.isin(positive + negative)]
        df['tf'] = df.index.isin(positive)
        df.sort_values(by=['avgFC'], inplace=True)
        return df

    @staticmethod
    def write_results(result_cfg, outdir):
        """

        :param result_cfg: results config file
        :param outdir: results output folder
        :return:
        """
        global RESULTS_FILE
        file_ext = '.tar.bz2'

        StaticMthods._create_tar(outdir + '/' + RESULTS_FILE + file_ext, outdir)
        generated_files = []
        for (dirpath, dirnames, filenames) in os.walk(outdir):
            generated_files.extend(dirnames)
            generated_files.extend(filenames)

        try:
            f = open(outdir + '/' + RESULTS_FILE + '.html', 'w')
            with open(result_cfg, 'r', encoding="utf-8") as cfgfile:
                cfg = json.load(cfgfile)
                rows = {'outdir': outdir, 'file_name': RESULTS_FILE + file_ext}
                f.write(''.join(cfg['header']).format(**rows))

                f.write(cfg['table_header'].format('I', 'pdf/Plotly images'))
                counter = 0
                for file_name, desc in cfg['images'].items():
                    if file_name in generated_files:
                        counter += 1
                        rows = {'count': counter, 'outdir': outdir, 'file_name': file_name,
                                'description': ' '.join(desc)}
                        f.write(cfg['table_row_images'].format(**rows))

                f.write(cfg['table_header'].format('II', 'Tab separated result files'))
                counter = 0
                for file_name, desc in cfg['files'].items():
                    if file_name in generated_files:
                        counter += 1
                        rows = {'count': counter, 'outdir': outdir, 'file_name': file_name,
                                'description': ' '.join(desc)}
                        f.write(cfg['table_row_files'].format(**rows))

                f.write(cfg['table_header'].format('III', 'Other algorithm output folders'))
                counter = 0
                for file_name, desc in cfg['folders'].items():
                    if file_name in generated_files:
                        StaticMthods._create_tar(outdir + '/' + file_name + file_ext, outdir + '/' + file_name)
                        counter += 1
                        rows = {'count': counter, 'outdir': outdir, 'file_name': file_name + file_ext,
                                'description': ' '.join(desc)}
                        f.write(cfg['table_row_folders'].format(**rows))

                f.write(cfg['footer'])
        except IOError as ioe:
            sys.exit('Can not create outfile:{}'.format(ioe.args[0]))

        return 0

    @staticmethod
    def _create_tar(outfile_name, source_dir):
        try:
            with tarfile.open(outfile_name, "w:bz2") as tar:
                tar.add(source_dir, arcname=os.path.basename(source_dir))
        except tarfile.TarError as te:
            log.info("error in tarfile generation:{}".format(te))
        except IOError as ioe:
            sys.exit('Error in tar file input folder file:{}'.format(ioe.args[0]))

        return 0
