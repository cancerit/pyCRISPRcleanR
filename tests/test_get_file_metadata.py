import archCompare.compareArchive as fc
import archCompare.staticMethods as sm
import pytest
import os

'''
written test to check codebase integrity
of archCompare
'''

class TestClass():
    testdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/')
    t_foldera=testdir+'testA'
    t_folderb=testdir+'testB'
    t_filea=testdir+'samplea.caveman_c.annot.vcf.gz'
    t_fileb=testdir+'testB.txt'
    t_tara=testdir+'testA.tar'
    t_tarb=testdir+'testB.tar'

    t_dirbamA=testdir+'testBamA'
    t_dirbamB=testdir+'testBamB'
    cwdpath=os.getcwd()
    configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
    t_json=configdir+'fileTypes.json'

    def test_get_file_metadata(self):
        my_dir_file=fc.ArchCompare(archive_a=self.t_foldera,archive_b=self.t_filea,json_config=self.t_json)
        assert('samplea.caveman_c.annot.vcf.gz','.vcf.gz') == my_dir_file._get_file_metadata(self.t_filea),'file metadata test ok'
        
if __name__ == '__main__':
  mytests=TestClass()
  mytests()
