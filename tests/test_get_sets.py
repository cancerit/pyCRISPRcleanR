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


    name_cmp_dict={'samplea.bam.bai': ['skipped', 'NoExtInJson'], 'samplea.caveman_c.annot.vcf.gz': ['compared', 'name'], 'samplea.bam': ['compared', 'name'], 'onlyinA.txt': ['skipped', 'onlyInA'], 'onlyinB.txt': ['skipped', 'onlyInB']}
    
    def test_get_sets(self):
        my_tar_file=fc.ArchCompare(archive_a=self.t_tara,archive_b=self.t_filea,json_config=self.t_json)
        assert self.name_cmp_dict == my_tar_file._get_sets_to_compare(self.format_dir_dictA,self.format_dir_dictB),'test_get_sets OK'
        
if __name__ == '__main__':
  mytests=TestClass()
  mytests()
