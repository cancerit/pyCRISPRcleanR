import logging.config
import os
import sys
from pyCRISPRcleanR.abstractCrispr import AbstractCrispr
from pyCRISPRcleanR.staticMethods import StaticMthods as SM
from pyCRISPRcleanR.plots import PlotData as PLT

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
logging.config.fileConfig(log_config)

log = logging.getLogger('cgpCRISPRcleanR')

'''
  This code run's Francesco's CRISPRcleanR algorithm implementation in python
'''


class CrisprCleanR(AbstractCrispr):
    """
       Main class , loads user defined paramters and files
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
        # check input files
        (input1, input2) = self.check_input()
        if input1 and input2:
            cldf = SM.combine_count_n_library(self.countfile, self.libfile, plot_flag=self.plot_data)
            cldf = SM.cleanup_data(cldf)
            cldf = SM.filter_data(cldf, controls, min_read_count)
            cldf, num_rep = SM.get_norm_count_n_fold_changes(cldf, controls, plot_flag=self.plot_data)
            cbs_dict = SM.run_cbs(cldf, cpus, sample)
            all_data = SM.process_segments(cbs_dict, ignored_genes, min_target_genes, controls, num_rep)

            if self.plot_data:
                cbs_dict_norm = SM.run_cbs(all_data, cpus, sample, fc_col="correctedFC")
                PLT.plot_segments(cbs_dict, cbs_dict_norm, sample)

            all_data.to_csv("corrected_counts_alldata.tsv", sep='\t', mode='w', header=True, index=True,
                            index_label='gRNA', doublequote=False)
            # print(all_data[all_data.gene == 'CCT8L2'])
        else:
            sys.exit('No valid input type provided')
