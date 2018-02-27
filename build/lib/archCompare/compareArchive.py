import json
import tarfile
import os
import sys
import tempfile
import re
import logging.config
from sys import stderr
from archCompare.abstractArchive import AbstractCompare
from archCompare.staticMethods import StaticMthods as sm

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
json_config = configdir + 'fileTypes.json'
logging.config.fileConfig(log_config)

log = logging.getLogger('compareArchive')

'''
Compare file, folder or archived data in .tar, .gz, .bz2 format
input any of the above
output -- list of items that are common in two compared archives
based on the MD5sum and/or data contents
'''


class ArchCompare(AbstractCompare):
    """
      Main class implements abstract class and
      its methods to check the inputype of a given file and loads parameters from user config
      file required in json format
    """

    def check_input(self):
        """ check input type of user supplied
          input filea and fileb keyword arguments
          and returns tupule of input file type
        """
        super().check_input()
        input_type = []
        for infile in (self.file_a, self.file_b):
            input_type.append(sm.input_checker(infile))
        return tuple(input_type)

    def get_config(self):
        """
          load parameters from json config file, if comparison type is not defined on command line
          uses default fron config file
        """
        super().get_config()
        try:
            if self.json_file is None:
                self.json_file = json_config
            with open(self.json_file, 'r') as cfgfile:
                self.cfg = json.load(cfgfile)
                self.prefix_ext = self.cfg['other_prm']['prefix_extension']
                if self.cmp_type is None:
                    self.cmp_type = ''.join(key for key, val in self.cfg['cmp_type'].items()
                                            if val.upper() == 'Y')
                    log.info('Using comaparison type from json config')
        except json.JSONDecodeError as jde:
            sys.exit('json error'.foramt(jde.args[0]))
        except FileNotFoundError as fne:
            sys.exit('Can not find json file:{}'.format(fne.args[0]))

    def _format_input(self, ftype, file_path):
        """
          accessory method to call other input formatting methods
        """
        if ftype == 'file':
            return self._format_file_input(file_path)
        elif ftype == 'tar':
            return self._format_tar_input(file_path)
        elif ftype == 'dir':
            return self._format_dir_input(file_path)
        else:
            log.erro('Undefined file format')
            return

    def _format_tar_input(self, file_path):
        """
          creates a diretory object of key = file name and values = [file paths, name, extensions]
          if requested comparison type is data then directory extracted to temp folder and
          input is then formatted
          by calling  _format_dir_input method
        """
        path_list = []
        list_for_prefix = []
        with tarfile.open(file_path, 'r') as tar:
            log.info('Processing tar file')
            if self.cmp_type in ['data', 'checksum']:
                log.info('Extracting tar file for data comparison, this might take a while ....')
                tmp_path = tempfile.mkdtemp(dir=".")
                tar.extractall(path=tmp_path)
                log.info(('Archive extraction completed at:', file_path))
                return self._format_dir_input(tmp_path)
            elif self.cmp_type == 'name':
                for tarinfo in tar:
                    if tarinfo.isreg():
                        name, ext = self._get_file_metadata(tarinfo.name)
                        path_list.append([tarinfo.name, name, ext])
                        if ext in self.prefix_ext:
                            list_for_prefix.append(name)
                prefix = os.path.commonprefix(list_for_prefix)
                return sm.process_list_to_dict(path_list, prefix)

    def _format_dir_input(self, file_path):
        """
          creates a diretory object of key = file name and values = [file paths, name, extensions]
        """
        path_list = []
        list_for_prefix = []
        log.info('Processing directory:{}'.format(file_path))
        for dirpath, _, files in os.walk(file_path):
            for filename in files:
                fullpath = os.path.join(dirpath, filename)
                name, ext = self._get_file_metadata(fullpath)
                path_list.append([fullpath, name, ext])
                if ext in self.prefix_ext:
                    list_for_prefix.append(name)
        prefix = os.path.commonprefix(list_for_prefix)
        return sm.process_list_to_dict(path_list, prefix)

    def _format_file_input(self, file_path):
        """
          creates a diretory object of key = file name and values = [file paths, name, extensions]
        """
        log.info('Processing file:{}'.format(file_path))
        name, ext = self._get_file_metadata(file_path)
        return sm.process_list_to_dict([[file_path, name, ext]], '')

    def _get_file_metadata(self, full_file_name):
        """
          takes file path as input and gives its path and processed extension
          # check second extension before .gz to determine file type [ e.g., .vcf.gz ]
        """
        (_, name) = os.path.split(full_file_name)
        (name_no_ext, first_ext) = os.path.splitext(name)
        if first_ext == '.gz':
            (_, second_ext) = os.path.splitext(name_no_ext)
            if second_ext in self.prefix_ext:
                first_ext = second_ext + first_ext
        return name, first_ext

    def _get_sets_to_compare(self, dictA, dictB):
        """
          peforms intersection and difference of diretory keys
          and outputs comparison of resulting sets as requested by user
          returns resuts dictionary containing filekey, comparison status and results if any
        """
        results_dict = {}
        common_files = list(set(dictA.keys()) & set(dictB.keys()))
        only_in_archiveA = list(set(dictA.keys()) - set(dictB.keys()))
        only_in_archiveB = list(set(dictB.keys()) - set(dictA.keys()))
        if self.cmp_type == 'name':
            for file_key in common_files:
                results_dict[file_key] = ['compared', 'name']
        else:
            results_dict = self._do_comparison(dictA, dictB, common_files)
        for file_key in only_in_archiveA:
            results_dict[file_key] = ['skipped', 'onlyInA']
        for file_key in only_in_archiveB:
            results_dict[file_key] = ['skipped', 'onlyInB']
        return results_dict

    def _do_comparison(self, dictA, dictB, common_files):
        """
          loops through dictionary and call explicit comaprsion methods as requested
          returns results dictionary containing filekey, comparison status and results if any
        """
        results_dict = {}
        json_data = self.cfg['extensions']
        checksum_type = self.cfg['other_prm']['checksum_type']
        for file_key in common_files:
            filea, _, ext_filea = (dictA[file_key])
            fileb, _, _ = (dictB[file_key])
            ext_dict = json_data.get(ext_filea, None)
            # retrieve json command for given extension
            if ext_dict is None:
                results_dict[file_key] = ['skipped', 'NoExtInJson']
            elif filea == fileb:
                log.info("Files have identical paths, skipping comaprison filea:{}fileb:{}".format(filea, fileb))
                results_dict[file_key] = ['skipped', 'IdenticalPath']
            elif self.cmp_type == 'checksum':
                log.info("performig checksum")
                result = sm.do_checksum_comaprison(checksum_type, filea=filea, fileb=fileb)
                results_dict[file_key] = ['compared', result]
            elif self.cmp_type == 'data':
                log.info("performig Data comparison")
                result = self._run_diff(ext_dict, ext_filea, filea=filea, fileb=fileb)
                results_dict[file_key] = ['compared', result]
            else:
                sys.exit('Unknown comparison type requested')
        return results_dict

    def _run_diff(self, ext_dict, ext, **kwargs):
        """ run comparison for given set of extension , additional methods could be added for different extension
          returns data [identical file content] or None [differences in file]
        """
        additional_prm = ext_dict.get('prm', None)
        cmd = ext_dict.get('cmd')
        exp_out = ext_dict.get('exp_out', None)
        log.info(("requested comparison type:", self.cmp_type))
        # add additional parametes
        if additional_prm is not None:
            for prm, val in additional_prm.items():
                kwargs[prm] = val
        (stdout, stderr) = sm.run_command(cmd.format(**kwargs))
        if stdout == 'Error':
            return 'Error'
        log.info(('ARGS:', kwargs, 'ERR:', stderr, 'OUT:', stdout))
        if re.search('^.vcf|^.vcf.gz', ext):
            return sm.get_vcf_diff(stdout, exp_out[0])
        else:
            if stderr:
                return
        return 'data'

    def run_comparison(self):
        """
          method to run the complete comparison
        """
        (typea, typeb) = self.check_input()
        dicta = self._format_input(typea, self.file_a)
        dictb = self._format_input(typeb, self.file_b)
        results = self._get_sets_to_compare(dicta, dictb)
        sm.format_results(results, dicta, dictb, self.outfile)
