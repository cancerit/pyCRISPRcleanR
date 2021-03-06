import logging
import json
import os
import sys
from pyCRISPRcleanR.abstractCrispr import AbstractCrispr
from pyCRISPRcleanR.staticMethods import StaticMthods as SM
from pyCRISPRcleanR.plots import PlotData as PLT
log = logging.getLogger(__name__)

'''
  This code run's Francesco's CRISPRcleanR algorithm implementation in python
'''


class CrisprCleanR(AbstractCrispr):
    """
        Main class , loads user defined parameters and files
        final data columns # sgRNA: guideRNA
        # <control sample count: raw 1..N replicates> : raw count
        # <treatment sample count: raw 1..N replicates> : raw count
        # gene: gene name as defined in the library file
        # chr: Chromosome name
        # start: gRNA start position
        # end: gRNA end position
        # <control sample count:normalised 1..N replicates> : Normalised count
        # <treatment sample count: normalised 1..N replicates> : Normalised count (postfixed _nc)
        # avgFC: average fold change values
        # BP: Base pair location ( used for DNAcopy analysis)
        # correction: correction factor
        # correctedFC: corrected foldchanges
        # <treatment sample count:corrected 1..N >: corrected count (postfixed _cc)

    """

    def check_input(self):
        """
           check input type and presence of user supplied
           input files
        """
        super().check_input()
        input_type = []
        for infile in (self.countfile, self.libfile):
            input_type.append(SM.input_checker(infile))
        return input_type

    def run_analysis(self):
        """
          method to run the analysis
        """
        controls = self.ncontrols
        min_read_count = self.minreads
        min_target_genes = self.mingenes
        ignored_genes = self.ignored_genes
        cpus = self.num_processors
        outdir = self.outdir
        gene_sig_dir = self.gene_sig_dir
        result_cfg = self.results_cfg
        iter = self.numiter

        global MAGECK_CMD
        if outdir:
            os.makedirs(outdir + '/mageckOut', exist_ok=True)
            os.makedirs(outdir + '/bagelOut', exist_ok=True)

        # check input files
        (input1, input2) = self.check_input()

        if input1 and input2:
            log.info("Running analysis, input file checks DONE.....")
            cldf = SM.combine_count_n_library(self.countfile, self.libfile, outdir=outdir)
            log.info("Count and library data combined.....")
            cldf = SM.filter_data(cldf, controls, min_read_count)
            log.info("Data filtering DONE.....")
            cldf, num_rep, norm_count_file, geneFC, sgRNAFC, fc, fcfile = \
                SM.get_norm_count_n_fold_changes(cldf, controls, outdir=outdir)
            log.info("Completed normalised count and fold change calculation .....")
            if self.run_mageck:
                norm_gene_summary = SM.run_mageck(norm_count_file, outfile_prefix=outdir + '/mageckOut/normCounts')
            if gene_sig_dir:
                ref_gene_list_dict = SM.load_signature_files(gene_sig_dir, cldf)
                if self.run_bagel:
                    log.info("Running Bagel on normalised fold changes .....")
                    cr = ','.join(map(lambda x: "%d" % x, range(1, fc.shape[1] - 1, 1)))
                    log.info("python BAGEL.py -i "+outdir+"/02_normalised_fold_changes.tsv -o "+outdir+"/bagelOut/normalised_FC_bagel.out -e "+gene_sig_dir+"/essential.txt -n "+gene_sig_dir+"/non_essential.txt -c "+str(cr)+" --numiter "+str(iter))
                    SM.run_bagel(fcfile, ref_gene_list_dict['essential_genes'],
                            ref_gene_list_dict['non_essential_genes'], cpus,
                            column_list=list(range(1, fc.shape[1] - 1, 1)), NUM_BOOTSTRAPS=iter,
                            outfilename=outdir + '/bagelOut/normalised_FC_bagel.out')
                    log.info("Completed Bagel on normalised fc .....")
                # ROC for sgRNA
                obs_pred_df = SM.get_obs_predictions(sgRNAFC, ref_gene_list_dict['essential_sgRNAs'],
                                                     ref_gene_list_dict['non_essential_sgRNAs'])
                PLT.roc_curve(obs_pred_df, data_type='sgRNA', saveto=outdir + '/04_roc_curve')
                PLT.pr_rc_curve(obs_pred_df, data_type='sgRNA', saveto=outdir + '/04_pr_rc_curve')
                # ROC for gene
                obs_pred_df = SM.get_obs_predictions(geneFC, ref_gene_list_dict['essential_genes'],
                                                     ref_gene_list_dict['non_essential_genes'])
                PLT.roc_curve(obs_pred_df, data_type='gene', saveto=outdir + '/05_roc_curve')
                PLT.pr_rc_curve(obs_pred_df, data_type='gene', saveto=outdir + '/05_pr_rc_curve')
                PLT.depletion_profile_with_gene_signature(geneFC, ref_gene_list_dict, obs_pred_df,
                                                          data_type='genes',
                                                          saveto=outdir + '/06_depletion_profile')
            # save normalised count and fold changes
            if self.runcrispr:
                cbs_dict = SM.run_cbs(cldf, cpus, fc_col='avgFC')
                log.info("CBS analysis completed  .....")
                all_data, corrected_count_file, crispr_fc, crispr_fc_file = SM.process_segments(cbs_dict,
                                                                                                ignored_genes,
                                                                                                min_target_genes,
                                                                                                controls, num_rep,
                                                                                                outdir=outdir)
                if gene_sig_dir:
                    essential, non_essential, other = SM.get_data_for_density_plot(all_data,
                                                                                   ref_gene_list_dict[
                                                                                       'essential_sgRNAs'],
                                                                                   ref_gene_list_dict[
                                                                                       'non_essential_sgRNAs'])
                    PLT.density_plot_ly(essential, non_essential, other,
                                        saveto=outdir + '/10_density_plots_pre_and_post_CRISPRcleanR')

                    if self.run_bagel:
                        log.info("Running Bagel on crisprcleanr corrected fold changes .....")
                        cr = ','.join(map(lambda x: "%d" % x, range(1, crispr_fc.shape[1] - 1, 1)))
                        log.info("python BAGEL.py -i "+outdir+"/04_crispr_cleanr_fold_changes.tsv -o "+outdir+"/bagelOut/CRISPRcleanR_FC_bagel.out -e "+gene_sig_dir+"/essential.txt -n "+gene_sig_dir+"/non_essential.txt -c "+str(cr)+" --numiter "+str(iter))
                        SM.run_bagel(crispr_fc_file, ref_gene_list_dict['essential_genes'],
                                ref_gene_list_dict['non_essential_genes'], cpus,
                                column_list=list(range(1, crispr_fc.shape[1] - 1, 1)), NUM_BOOTSTRAPS=iter,
                                outfilename=outdir + '/bagelOut/CRISPRcleanR_FC_bagel.out')
                        log.info("Completed Bagel on crisprcleanr corrected fold changes .....")

                log.info("Processed CBS segments  .....")
                SM._print_df(all_data, outdir + "/05_alldata.tsv")
                cbs_dict_norm = SM.run_cbs(all_data, cpus, fc_col="correctedFC")
                log.info("CBS analysis on normalised fold changes completed.....")
                PLT.plot_segments(cbs_dict, cbs_dict_norm, outdir=outdir)
                log.info("Done plots.....")

                if self.run_mageck:
                    log.info("Running Mageck.....")
                    corrected_gene_summary = SM.run_mageck(corrected_count_file,
                                                           outfile_prefix=outdir + '/mageckOut/correctedCounts')

                    PLT.impact_on_phenotype(norm_gene_summary, corrected_gene_summary,
                                            saveto=outdir + '/11_impact_on_phenotype')

            SM.write_results(result_cfg, outdir)
            log.info("Analysis completed successfully.....")
        else:
            sys.exit('Input data is not in required format, see inputFormat in README file')
