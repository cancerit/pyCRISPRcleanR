import pyCRISPRcleanR.formatInput as fc
import pyCRISPRcleanR.staticMethods as sm
import pytest
import os

'''
written test to check codebase integrity
of archCompare
'''

class TestClass():
    testdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/')
    t_countfile=testdir + 'KY_libv1.0_chr22.tsv.gz'
    t_libfile = testdir + 'HT-29_counts_chr22.tsv.gz'
    t_sig_dir = testdir + 'ref_genes.tar.gz'

    cwdpath=os.getcwd()

    def test_file_input(self):
        # check input type function
        my_tar_tar=fc.CrisprCleanR(countfile=self.t_countfile, libfile=self.t_libfile, gene_signatures=self.t_sig_dir)
        assert ['file','file'] == my_tar_tar.check_input(),'tar n tar test OK'


if __name__ == '__main__':
  mytests=TestClass()
  mytests()
