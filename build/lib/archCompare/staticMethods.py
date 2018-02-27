import sys
import os
import tarfile
from subprocess import Popen, PIPE, STDOUT
import shlex
import re
import logging.config

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
json_config = configdir + 'fileTypes.json'
logging.config.fileConfig(log_config)

log = logging.getLogger('compareArchive')


class StaticMthods(object):
    """ Static methosds for sommon tasks """

    def __init__(self):
        super().__init__()

    def input_checker(infile):
        """
          checks user input file and returns it's type
        """
        try:
            if tarfile.is_tarfile(infile):
                log.info(("input is an archive:", infile))
                return 'tar'
            else:
                log.info(("input is a file:", infile))
                return 'file'
        except IsADirectoryError:
            return 'dir'
        except IOError as ioe:
            sys.exit('Error in reading input file:{}'.format(ioe.args[0]))

    def process_list_to_dict(path_list, prefix):
        """
          creates dictionary of paths [relative path] in an archive
          where filename prefix replaced to create keys that are identical in two dictionaries
          path_list = [file_path, name, ext]
        """
        path_dict = {}
        for meta_list in path_list:
            meta_list[0] = os.path.realpath(meta_list[0])
            path_dict[meta_list[1].replace(prefix, '', 1)] = meta_list
        return path_dict

    def run_command(cmd):
        """ runs command in a shell, returns stdout and exit code"""
        args = shlex.split(cmd)
        try:
            cmd_obj = Popen(args, stdin=None, stdout=PIPE, stderr=STDOUT,
                            shell=False, universal_newlines=True, bufsize=1)
            (stdout, stderr) = cmd_obj.communicate()
            return stdout, cmd_obj.returncode
        except OSError as oe:
            log.error(("Unable to run command", cmd, oe.args))
            return 'Error', 'Error'

    def format_results(results, dicta, dictb, outfile):
        """
          formats output results as tab separated files
        """
        if outfile is not None:
            f = open('results.tsv', 'w')
            f.write('Filea\tFileb\tStatus\tSimilarityBy\n')

        print('Filea\tFileb\tStatus\tSimilarityBy\n')
        for key, value in results.items():
            if value[1] is None:
                value[1] = 'differ'
            if outfile is None:
                print("{}\t{}\t{}\t{}\n".format(dicta.get(key, ['NA'])[0],
                                                dictb.get(key, ['NA'])[0], value[0], value[1]))
            else:
                f.write("{}\t{}\t{}\t{}\n".format(dicta.get(key, ['NA'])[0],
                                                  dictb.get(key, ['NA'])[0], value[0], value[1]))

    def do_checksum_comaprison(prog, **kwargs):
        """
          peforms checksum comparison and
          returns checksum [identical] or None [different checksum values]
        """
        cmd = r'{checksum} {filea} {fileb}'
        kwargs['checksum'] = prog
        (stdout, stderr) = StaticMthods.run_command(cmd.format(**kwargs))
        out_msg = re.split('\n|\s', stdout)
        if stdout == 'Error':
            return 'Error'
        elif out_msg[0] == out_msg[3]:
            return 'checksum'
        else:
            return

    def get_vcf_diff(stdout, exp_out):
        """
          compres two vcf file and returns
          data [identical file content] or None [differences in file]
        """
        out_msg = stdout.split("\n")
        for line in out_msg:
            if not line:  # skip empty lines
                continue
            try:
                match = re.search(exp_out, line)
                if match is None:  # if match not found return true
                    return 'data'
                elif match[0] == '##contig=':  # matches only when contigs are different in two vcf files
                    return
            except re.error:
                print('Error in regular expression:{}'.format(re.error))
        return 'data'
