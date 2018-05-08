import logging
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
        # correctedFC: corrected foldchange values
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
        sample = self.sample
        cpus = self.num_processors
        outdir = self.outdir
        if outdir:
            os.makedirs(outdir, exist_ok=True)

        # check input files
        (input1, input2) = self.check_input()

        if input1 and input2:
            log.info("Running analysis, input file checks DONE.....")
            cldf = SM.combine_count_n_library(self.countfile, self.libfile,
                                              plot_flag=self.plot_data, outdir=outdir)
            log.info("Count and library data combined.....")
            cldf = SM.filter_data(cldf, controls, min_read_count)
            log.info("Data filtering DONE.....")
            cldf, num_rep = SM.get_norm_count_n_fold_changes(cldf, controls,
                                                             plot_flag=self.plot_data, outdir=outdir)
            log.info("Completed normalised count and fold change calculation .....")

            # save normalised count and fold changes
            if self.runcrispr:
                cbs_dict = SM.run_cbs(cldf, cpus, sample)
                log.info("CBS analysis completed  .....")
                all_data = SM.process_segments(cbs_dict, ignored_genes, min_target_genes, controls, num_rep,
                                               outdir=outdir)
                log.info("Processed CBS segments  .....")
                SM._print_df(all_data, outdir + "/alldata.tsv")
                if self.plot_data:
                    cbs_dict_norm = SM.run_cbs(all_data, cpus, sample, fc_col="correctedFC")
                    log.info("CBS analysis on normalised fold changes completed.....")
                    PLT.plot_segments(cbs_dict, cbs_dict_norm, sample, outdir=outdir)
                    log.info("Done plots.....")
            log.info("Analysis completed successfully.....")
            print("Analysis completed successfully.....")
            # print(all_data[all_data.gene == 'CCT8L2'])
        else:
            sys.exit('Input data is not in required format, see inputFormat in README file')
