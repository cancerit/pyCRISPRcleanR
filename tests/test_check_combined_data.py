import pyCRISPRcleanR.formatInput as fc
import pyCRISPRcleanR.staticMethods as sm
import pyCRISPRcleanR.plots as PLT
import pandas as pd
import pytest
import os, tempfile
import filecmp

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
        expected_norm_gene_summary = testdir + 'mageck_gene_summary_norm.txt'
        expected_corrected_gene_summary = testdir + 'mageck_gene_summary_corrected.txt'
        controls = 1
        min_read_count = 30
        min_target_genes = 3
        ignored_genes = []
        sample = 'mytest'
        gene_sig_dir= testdir + '/ref_genes/'
        cpus = 1
        outdir=tempfile.mkdtemp(dir=".")
        if outdir:
            os.makedirs(outdir + '/mageckOut', exist_ok=True)
        # check input type function
        mystatic_obj = sm.StaticMthods()
        myPLT = PLT.PlotData()
        cldf=mystatic_obj.combine_count_n_library(t_countfile, t_libfile, outdir=outdir)
        assert (2072, 8) == cldf.shape,'combined_count_n_lib'
        cldf=mystatic_obj.filter_data(cldf, controls, min_read_count)
        assert (2038, 8) == cldf.shape, 'filter_data'
        cldf,no_rep, norm_counts_file, geneFC, sgRNAFC = mystatic_obj.get_norm_count_n_fold_changes(cldf,controls,outdir=outdir)
        assert (2038, 14) == cldf.shape, 'get_norm_count_n_fold_changes'
        assert 3 == no_rep, 'number of replicates'
        ref_gene_list_dict = mystatic_obj.load_signature_files(gene_sig_dir, cldf)

        obs_pred_df = mystatic_obj.get_obs_predictions(sgRNAFC, ref_gene_list_dict['essential_sgRNAs'],
                                             ref_gene_list_dict['non_essential_sgRNAs'])
        assert (94,2) == obs_pred_df.shape, 'get data for roc'
        myPLT.roc_curve(obs_pred_df, data_type='sgRNA', saveto=outdir + '/roc_curve')
        myPLT.pr_rc_curve(obs_pred_df, data_type='sgRNA', saveto=outdir + '/pr_rc_curve')
        cbs_dict=mystatic_obj.run_cbs(cldf, cpus, sample)
        assert "dict_keys([22])" == "{}".format(cbs_dict.keys()), 'run_cbs'
        alldata, corrected_counts_file = mystatic_obj.process_segments(cbs_dict, ignored_genes, min_target_genes, controls, no_rep,outdir=outdir)
        # only used to create final test results
        #alldata.to_pickle('pickled_df_HT-29.pkl', compression='gzip', protocol=-1)
        expected_df=pd.read_pickle(picke_file, compression='gzip')
        result=expected_df.equals(alldata)
        assert (2038, 24) == alldata.shape, 'process_segments'
        # test code no checks
        norm_gene_summary, corrected_gene_summary = mystatic_obj.run_mageck(norm_counts_file,corrected_counts_file,outdir=outdir)
        assert filecmp.cmp(expected_norm_gene_summary,norm_gene_summary,shallow=True) , 'mageck files are identical'
        assert filecmp.cmp(expected_corrected_gene_summary, corrected_gene_summary, shallow=True), 'mageck files are identical'
        cbs_dict_norm = mystatic_obj.run_cbs(alldata, cpus, sample, fc_col="correctedFC")
        myPLT.plot_segments(cbs_dict, cbs_dict_norm, sample, outdir=outdir)
        #assert True == result, 'process_segments: check results'

if __name__ == '__main__':
  mytests=TestClass()
  mytests()
