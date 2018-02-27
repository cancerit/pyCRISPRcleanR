from archCompare.compareArchive import ArchCompare
import sys
import os
import argparse
import pkg_resources

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
version = pkg_resources.require("archCompare")[0].version


def main():
    usage = "\n %prog [options] -a archive_a.tar -b archive_b.tar -j fileTypes.json -c name"

    optParser = argparse.ArgumentParser(prog='cgpCompare')
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-a", "--archive_a", type=str, dest="archive_a", required=True,
                          default="", help="archive path for a folder, file or a tar data")

    required.add_argument("-b", "--archive_b", type=str, dest="archive_b", required=True,
                          default="", help="archive path for a folder, file or a tar data")

    optional.add_argument("-j", "--json_config", type=str, dest="json_config",
                          default=None, help="path to json config file")

    optional.add_argument("-o", "--outfile", type=str, dest="outfile",
                          default=None, help="path to outfile file, STOUT if not provided")

    optional.add_argument("-c", "--cmp_type", type=str, dest="cmp_type", required=False,
                          default=None, help="Compariosn type to perform [ \
                           name: compares archives using file names, \
                           checksum: perform comparsion based on checksum tool defined in json config file \
                           data: does full data comparison based on tools defined for each extension type \
                           Note- command line option if set overrides default value in json config file ]")

    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + version)
    optional.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True)

    optParser._action_groups.append(optional)

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)
    opts = optParser.parse_args()
    if not (opts.archive_a or opts.archive_b):
        sys.exit('\nERROR Arguments required\n\tPlease run: cgpCompare --help\n')
    mycomp = ArchCompare(archive_a=opts.archive_a,
                         archive_b=opts.archive_b,
                         json_config=opts.json_config,
                         cmp_type=opts.cmp_type,
                         outfile=opts.outfile)
    mycomp.run_comparison()
