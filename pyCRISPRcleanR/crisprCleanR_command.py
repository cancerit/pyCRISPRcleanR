from pyCRISPRcleanR.formatInput import CrisprCleanR
import sys
import os
import argparse
import pkg_resources
import logging.config

configdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config/')
results_json = configdir + 'results.json'
log_config = configdir + 'logging.conf'
ref_genes = configdir + 'ref_genes.tar.gz'
logging.config.fileConfig(log_config)
log = logging.getLogger(__name__)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
version = pkg_resources.require("pyCRISPRcleanR")[0].version


def main():  # pragma: no cover
    usage = "\n %prog [options] -f counts.tsv -l library.tsv"

    optParser = argparse.ArgumentParser(prog='pyCRISPRcleanR',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = optParser._action_groups.pop()
    required = optParser.add_argument_group('required arguments')

    required.add_argument("-f", "--countfile", type=str, dest="countfile", required=True,
                          default="", help="sgRNA raw count file, accepts compressed file")

    required.add_argument("-l", "--libfile", type=str, dest="libfile", required=True,
                          default="", help="sgRNA library file, accepts compressed file")

    optional.add_argument("-mr", "--minreads", type=int, dest="minreads", required=False,
                          default=30, help="minimum read count in control sample \
                          to be used for filtering ")

    optional.add_argument("-mg", "--mingenes", type=int, dest="mingenes", required=False,
                          default=3, help="minimum number of genes in a CNV segment to\
                          consider it for count normalization ")

    optional.add_argument("-ig", "--ignored_genes", nargs='*', type=str, dest="ignored_genes", required=False,
                          default=[], help="space separate list of ignored genes")

    optional.add_argument("-nc", "--ncontrols", type=int, dest="ncontrols", required=False,
                          default=1, help="Number of control samples in raw count file [ \
                           Note: at least one control sample is required ]")

    optional.add_argument("-np", "--num_processors", type=int, dest="num_processors", required=False,
                          default=1, help="Number of processors to use for parallel jobs")

    optional.add_argument("-cc", "--crispr_cleanr", action='store_true', dest="crispr_cleanr",
                          help="flag to run CRISPRcleanR")

    optional.add_argument("-gs", "--gene_signatures", type=str, dest="gene_signatures", required=False,
                          default=ref_genes, help="Path to a directory or tar archive \
                          containing files for signature genes")

    optional.add_argument("-mk", "--run_mageck", action='store_true', dest="run_mageck",
                          help="flag to run MAGeCK")

    optional.add_argument("-bl", "--run_bagel", action='store_true', dest="run_bagel",
                          help="flag to run BAGEL")

    optional.add_argument("-N", "--numiter", type=int, dest="numiter", required=False,
                          default=1000, help="Number of bootstrap iterations for BAGEL (default 1000)")

    optional.add_argument("-o", "--outdir", type=str, dest="outdir",
                          default='./', help="path to output folder ")

    optional.add_argument("-results", "--results_cfg", type=str, dest="results_cfg",
                          default=results_json, help=argparse.SUPPRESS)

    optional.add_argument("-v", "--version", action='version', version='%(prog)s ' + version)
    optional.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True)

    optParser._action_groups.append(optional)

    if len(sys.argv) == 1:
        optParser.print_help()
        log.debug("Missing one or more required arguments in the command, exiting......")
        sys.exit("Missing one or more required arguments in the command, exiting......")
    opts = optParser.parse_args()
    if not (opts.countfile or opts.libfile):
        log.debug('ERROR Arguments required \n Please run: pyCRISPRcleanR --help')
        sys.exit('\nERROR Arguments required\n\tPlease run: pyCRISPRcleanR --help\n')
    log.debug('Analysis started....')
    log.info(opts)
    mycrispr = CrisprCleanR(**vars(opts))
    mycrispr.run_analysis()


if __name__ == '__main__':
    main()
