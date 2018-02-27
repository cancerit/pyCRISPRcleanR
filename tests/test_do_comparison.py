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
    cmp_type='data'
    file_dict={'samplea.caveman_c.annot.vcf.gz':[t_filea,'samplea.caveman_c.annot.vcf.gz','.vcf.gz']}

    format_dir_bamdiffA={'samplea.bam': [t_dirbamA+'/samplea.bam', 'samplea.bam', '.bam']}
    format_dir_bamdiffB={'samplea.bam': [t_dirbamB+'/samplea.bam', 'samplea.bam', '.bam']}
    common_files_bamdiff=['samplea.bam']
    bamdiff_result={'samplea.bam': ['compared', None]}


    format_tar_dictA={
                      'bam': [cwdpath+'/testA/samplea.bam', 'samplea.bam', '.bam'],
                      'bam.bai': [cwdpath+'/testA/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinA.txt': [cwdpath+'/testA/onlyinA.txt', 'onlyinA.txt', '.txt'],
                      'caveman_c.annot.vcf.gz': [cwdpath+'/testA/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz']}


    format_tar_dictB={
                      'bam': [cwdpath+'/testB/samplea.bam', 'samplea.bam', '.bam'],
                      'caveman_c.annot.vcf.gz': [cwdpath+'/testB/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz'],
                      'bam.bai': [cwdpath+'/testB/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinB.txt': [cwdpath+'/testB/onlyinB.txt', 'onlyinB.txt', '.txt']}



    format_dir_dictA={
                      'bam': [t_foldera+'/samplea.bam', 'samplea.bam', '.bam'],
                      'bam.bai': [t_foldera+'/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinA.txt': [t_foldera+'/onlyinA.txt', 'onlyinA.txt', '.txt'],
                      'caveman_c.annot.vcf.gz': [t_foldera+'/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz']}


    format_dir_dictB={
                      'bam': [t_folderb+'/samplea.bam', 'samplea.bam', '.bam'],
                      'caveman_c.annot.vcf.gz': [t_folderb+'/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz'],
                      'bam.bai': [t_folderb+'/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinB.txt': [t_folderb+'/onlyinB.txt', 'onlyinB.txt', '.txt']}

    common_inAB=['samplea.bam', 'samplea.bam.bai', 'samplea.caveman_c.annot.vcf.gz']
    only_inA=['onlyinA.txt']
    only_inB=['onlyinB.txt']

    common_files=['bam', 'bam.bai', 'caveman_c.annot.vcf.gz']

    name_cmp_dict={'bam.bai': ['compared', 'name'], 'caveman_c.annot.vcf.gz': ['compared', 'name'], 'bam': ['compared', 'name'], 'onlyinA.txt': ['skipped', 'onlyInA'], 'onlyinB.txt': ['skipped', 'onlyInB']}

    data_cmp_dict={'bam': ['compared', 'data'], 'bam.bai': ['skipped', 'NoExtInJson'], 'caveman_c.annot.vcf.gz': ['compared', None]}
    checksum_cmp_dict={'bam': ['compared', 'checksum'], 'bam.bai': ['skipped', 'NoExtInJson'], 'caveman_c.annot.vcf.gz': ['compared', None]}

    @pytest.mark.skipif("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true", reason="Skipping this test on Travis CI.")
    def test_do_comparison(self):
        self.maxDiff = None
        my_tar_tar_cmp=fc.ArchCompare(archive_a=self.t_tara,archive_b=self.t_tarb,json_config=self.t_json,cmp_type='data')
        assert self.data_cmp_dict == my_tar_tar_cmp._do_comparison(self.format_dir_dictA,self.format_dir_dictB,self.common_files),'test_do_comparison OK'

    @pytest.mark.skipif("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true", reason="Skipping this test on Travis CI.")
    def test_bamdiff(self):
        self.maxDiff = None
        my_dir_dir_bam=fc.ArchCompare(archive_a=self.t_dirbamA,archive_b=self.t_dirbamB,json_config=self.t_json,cmp_type='data')
        assert self.bamdiff_result == my_dir_dir_bam._do_comparison(self.format_dir_bamdiffA,self.format_dir_bamdiffB,self.common_files_bamdiff),'test_do_BamComparison OK'

    def test_checksum(self):
        self.maxDiff = None
        my_dir_dir_checksum=fc.ArchCompare(archive_a=self.t_foldera,archive_b=self.t_folderb,json_config=self.t_json,cmp_type='checksum')
        assert self.checksum_cmp_dict == my_dir_dir_checksum._do_comparison(self.format_dir_dictA,self.format_dir_dictB,self.common_files),'test_do_BamComparison OK'

if __name__ == '__main__':
  mytests=TestClass()
  mytests()
