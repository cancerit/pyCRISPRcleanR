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
                      'samplea.bam': [cwdpath+'/testA/samplea.bam', 'samplea.bam', '.bam'],
                      'samplea.bam.bai': [cwdpath+'/testA/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinA.txt': [cwdpath+'/testA/onlyinA.txt', 'onlyinA.txt', '.txt'],
                      'samplea.caveman_c.annot.vcf.gz': [cwdpath+'/testA/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz']}


    format_tar_dictB={
                      'samplea.bam': [cwdpath+'/testB/samplea.bam', 'samplea.bam', '.bam'],
                      'samplea.caveman_c.annot.vcf.gz': [cwdpath+'/testB/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz'],
                      'samplea.bam.bai': [cwdpath+'/testB/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinB.txt': [cwdpath+'/testB/onlyinB.txt', 'onlyinB.txt', '.txt']}



    format_dir_dictA={
                      'samplea.bam': [t_foldera+'/samplea.bam', 'samplea.bam', '.bam'],
                      'samplea.bam.bai': [t_foldera+'/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinA.txt': [t_foldera+'/onlyinA.txt', 'onlyinA.txt', '.txt'],
                      'samplea.caveman_c.annot.vcf.gz': [t_foldera+'/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz']}


    format_dir_dictB={
                      'samplea.bam': [t_folderb+'/samplea.bam', 'samplea.bam', '.bam'],
                      'samplea.caveman_c.annot.vcf.gz': [t_folderb+'/vcf_data/samplea.caveman_c.annot.vcf.gz', 'samplea.caveman_c.annot.vcf.gz', '.vcf.gz'],
                      'samplea.bam.bai': [t_folderb+'/samplea.bam.bai', 'samplea.bam.bai', '.bai'],
                      'onlyinB.txt': [t_folderb+'/onlyinB.txt', 'onlyinB.txt', '.txt']}

    common_inAB=['samplea.bam', 'samplea.bam.bai', 'samplea.caveman_c.annot.vcf.gz']
    only_inA=['onlyinA.txt']
    only_inB=['onlyinB.txt']

    common_files=['bam', 'bam.bai', 'caveman_c.annot.vcf.gz']

    name_cmp_dict={'bam.bai': ['compared', 'name'], 'caveman_c.annot.vcf.gz': ['compared', 'name'], 'bam': ['compared', 'name'], 'onlyinA.txt': ['skipped', 'onlyInA'], 'onlyinB.txt': ['skipped', 'onlyInB']}

    data_cmp_dict={'bam': ['compared', 'data'], 'bam.bai': ['skipped', 'NoExtInJson'], 'caveman_c.annot.vcf.gz': ['compared', None]}
    checksum_cmp_dict={'bam': ['compared', 'checksum'], 'bam.bai': ['skipped', 'NoExtInJson'], 'caveman_c.annot.vcf.gz': ['compared', None]}


    def test_format_file_input(self):
        my_dir_file=fc.ArchCompare(archive_a=self.t_foldera,archive_b=self.t_filea,json_config=self.t_json,cmp_type=self.cmp_type)
        assert self.file_dict == my_dir_file._format_file_input(self.t_filea),'test_format_file_input OK'

    def test_format_dir_input(self):
        self.maxDiff = None
        my_dir_file=fc.ArchCompare(archive_a=self.t_foldera,archive_b=self.t_filea,json_config=self.t_json,cmp_type=self.cmp_type)
        assert self.format_dir_dictA == my_dir_file._format_dir_input(self.t_foldera),'test_format_dir_input OK'

    def test_format_tar_input(self):
        self.maxDiff = None
        my_tar_tar=fc.ArchCompare(archive_a=self.t_tara,archive_b=self.t_tarb,json_config=self.t_json,cmp_type='name')
        assert self.format_tar_dictA == my_tar_tar._format_tar_input(self.t_tara),'test_format_tar_inputA OK'
        assert self.format_tar_dictB == my_tar_tar._format_tar_input(self.t_tarb),'test_format_tar_inputB OK'

if __name__ == '__main__':
  mytests=TestClass()
  mytests()
