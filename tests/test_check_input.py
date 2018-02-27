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
   
    def test_dir_file_type(self):
        # check input type function
        my_dir_file=fc.ArchCompare(archive_a=self.t_foldera,archive_b=self.t_filea,json_config=self.t_json)
        assert my_dir_file.check_input() == ('dir','file') , 'Directory n file test OK'
        my_dir_tar=fc.ArchCompare(archive_a=self.t_foldera,archive_b=self.t_tara,json_config=self.t_json)
        assert ('dir','tar')== my_dir_tar.check_input(),'directory n tar test OK'
        my_tar_tar=fc.ArchCompare(archive_a=self.t_tarb,archive_b=self.t_tara,json_config=self.t_json)
        assert ('tar','tar') == my_tar_tar.check_input(),'tar n tar test OK'
        my_dir_dir=fc.ArchCompare(archive_a=self.t_foldera,archive_b=self.t_folderb,json_config=self.t_json)
        assert ('dir','dir')== my_dir_dir.check_input(),'directory n directory test OK'
        my_file_file=fc.ArchCompare(archive_a=self.t_filea,archive_b=self.t_fileb,json_config=self.t_json)
        assert ('file','file') == my_file_file.check_input(),'file n file test OK'

if __name__ == '__main__':
  mytests=TestClass()
  mytests()
