from cgpCRISPRcleanR.formatInput import CrisprCleanR
import sys
import os
import argparse
import pkg_resources

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
version = pkg_resources.require("cgpCrisprCleanR")[0].version


def main():
    usage = "\n %prog [options] -f counts.tsv -l library.tsv"

    optParser = argparse.ArgumentParser(prog='runCrisprCleanR')
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-f", "--countfile", type=str, dest="countfile", required=True,
                          default="", help="sgRNA raw count file")

    required.add_argument("-l", "--libfile", type=str, dest="libfile", required=True,
                          default="", help="sgRNA library file")

    optional.add_argument("-e", "--expname", type=str, dest="expname", required=False,
                          default='crisprExperiment', help="name of the experiment")

    optional.add_argument("-r", "--minreads", type=str, dest="minreads", required=False,
                          default=30, help="minimum read count in control sample \
                          to be used for filtering ")

    optional.add_argument("-g", "--mingenes", type=str, dest="mingenes", required=False,
                          default=3, help="minimum number of genes in a CNV segment to\
                          consider it for count normalization ")

    optional.add_argument("-c", "--ncontrols", type=str, dest="ncontrols", required=False,
                          default=1, help="Number of control samples in raw count file [ \
                           Note: at least one control sample is required ]")

    optional.add_argument("-s", "--sample", type=str, dest="sample", required=False,
                          default='mysample', help="sample name in counts file")

    optional.add_argument("-o", "--outdir", type=str, dest="outdir",
                          default=None, help="path to output folder")

    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + version)
    optional.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True)

    optParser._action_groups.append(optional)

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)
    opts = optParser.parse_args()
    if not (opts.countfile or opts.libfile):
        sys.exit('\nERROR Arguments required\n\tPlease run: cgpCrispCleanR --help\n')
    mycrispr = CrisprCleanR(**vars(opts))
    mycrispr.run_analysis()
