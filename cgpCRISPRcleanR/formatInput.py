import logging.config
import os
from cgpCRISPRcleanR.abstractCrispr import AbstractCrispr
from cgpCRISPRcleanR.staticMethods import StaticMthods as SM

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
json_config = configdir + 'fileTypes.json'
logging.config.fileConfig(log_config)

log = logging.getLogger('cgpCRISPRcleanR')

'''
   run's Francesco's CRISPRcleanR algorithm implemented in pyth    on
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
            print("File:{} is type:{}".format(infile, input_type))
        return None

    def run_analysis(self):
        """
          method to run the analysis
        """
        controls = self.ncontrols
        min_read_count = self.minreads
        min_target_genes = self.mingenes
        ignored_genes = []  # should come from command line ???
        sample = 'mysample'
        cpus = 4
        # check input files
        self.check_input()
        cldf = SM.combine_count_n_library(self.countfile, self.libfile)
        cldf = SM.cleanup_data(cldf)
        cldf = SM.filter_data(cldf, controls, min_read_count)
        (cldf, num_rep) = SM.get_norm_count_n_fold_changes(cldf, controls)
        cbs_dict = SM.run_cbs(cldf, cpus, sample)
        all_data = SM.process_segments(cbs_dict, ignored_genes, min_target_genes, controls, num_rep)
        all_data.to_csv("corrected_counts_alldf_v9nc.tsv", sep='\t', mode='w', header=True, index=True,
                        index_label='gRNA', doublequote=False)
        print(all_data[all_data.gene == 'A1BG'])
