import cgpCRISPRcleanR.formatInput as fc
import cgpCRISPRcleanR.staticMethods as sm
import pandas as pd
import pytest
import os

'''
written test to check codebase integrity
of archCompare
'''



class TestClass():
    pass

    def test_static_methods(self):
        testdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/')
        t_libfile = testdir + 'KY_libv1.0_chr22.tsv.gz'
        t_countfile = testdir + 'HT-29_counts_chr22.tsv.gz'
        picke_file = testdir + 'pickled_df_HT-29.pkl'
        controls = 1
        min_read_count = 30
        min_target_genes = 3
        ignored_genes = []
        sample = 'mytest'
        cpus = 1
        # check input type function
        mystatic_obj=sm.StaticMthods()
        cldf=mystatic_obj.combine_count_n_library(t_countfile, t_libfile)
        assert (2072, 12) == cldf.shape,'combined_count_n_lib'
        cldf=mystatic_obj.cleanup_data(cldf)
        assert (2072, 8) == cldf.shape, 'cleanup_data'
        cldf=mystatic_obj.filter_data(cldf, controls, min_read_count)
        assert (2038, 8) == cldf.shape, 'filter_data'
        cldf,no_rep=mystatic_obj.get_norm_count_n_fold_changes(cldf,controls)
        assert (2038, 14) == cldf.shape, 'get_norm_count_n_fold_changes'
        assert 3 == no_rep, 'number of replicates'
        cbs_dict=mystatic_obj.run_cbs(cldf, cpus, sample)
        assert "dict_keys([22])" == "{}".format(cbs_dict.keys()), 'run_cbs'
        alldata = mystatic_obj.process_segments(cbs_dict, ignored_genes, min_target_genes, controls, no_rep)
        # only used to create final test results
        #alldata.to_pickle('pickled_df_HT-29.pkl', compression='gzip', protocol=-1)
        expected_df=pd.read_pickle(picke_file, compression='gzip')
        result=expected_df.equals(alldata)
        assert (2038, 19) == alldata.shape, 'process_segments'
        assert True == result, 'process_segments: check results'

if __name__ == '__main__':
  mytests=TestClass()
  mytests()
