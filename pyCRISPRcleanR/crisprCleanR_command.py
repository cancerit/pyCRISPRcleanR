from pyCRISPRcleanR.formatInput import CrisprCleanR
import sys
import os
import argparse
import pkg_resources
import logging.config

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
log_config = configdir + 'logging.conf'
logging.config.fileConfig(log_config)
log = logging.getLogger(__name__)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
version = pkg_resources.require("pyCRISPRcleanR")[0].version


def main():
    usage = "\n %prog [options] -f counts.tsv -l library.tsv"

    optParser = argparse.ArgumentParser(prog='pyCRISPRCleanR',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-f", "--countfile", type=str, dest="countfile", required=True,
                          default="", help="sgRNA raw count file, accepts compressed file")

    required.add_argument("-l", "--libfile", type=str, dest="libfile", required=True,
                          default="", help="sgRNA library file, accepts compressed file")

    optional.add_argument("-e", "--expname", type=str, dest="expname", required=False,
                          default='myexperiment', help="name of the experiment")

    optional.add_argument("-mr", "--minreads", type=str, dest="minreads", required=False,
                          default=30, help="minimum read count in control sample \
                          to be used for filtering ")

    optional.add_argument("-mg", "--mingenes", type=str, dest="mingenes", required=False,
                          default=3, help="minimum number of genes in a CNV segment to\
                          consider it for count normalization ")

    optional.add_argument("-ig", "--ignored_genes", nargs='*', type=str, dest="ignored_genes", required=False,
                          default=[], help="space separate list of ignored genes")

    optional.add_argument("-nc", "--ncontrols", type=str, dest="ncontrols", required=False,
                          default=1, help="Number of control samples in raw count file [ \
                           Note: at least one control sample is required ]")

    optional.add_argument("-s", "--sample", type=str, dest="sample", required=False,
                          default='mysample', help="sample name in counts file")

    optional.add_argument("-np", "--num_processors", type=int, dest="num_processors", required=False,
                          default=1, help="Number of processors to use for parallel jobs")

    optional.add_argument("-pl", "--plot_data", type=str, dest="plot_data", required=False,
                          default=None, help="Generate pdf and interactive plotly images [y or 1]")

    optional.add_argument("-o", "--outdir", type=str, dest="outdir",
                          default='./', help="path to output folder ")

    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + version)
    optional.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True)

    optParser._action_groups.append(optional)

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)
    opts = optParser.parse_args()
    if not (opts.countfile or opts.libfile):
        sys.exit('\nERROR Arguments required\n\tPlease run: pyCRISPRCleanR --help\n')
        log.debug('ERROR Arguments required \n Please run: pyCRISPRCleanR --help')
    log.debug('Analysis started....')

    mycrispr = CrisprCleanR(**vars(opts))
    mycrispr.run_analysis()


if __name__ == '__main__':
    main()
