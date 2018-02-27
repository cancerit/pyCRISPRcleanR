import sys,os
import argparse
import pkg_resources 

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

from  archCompare.compareArchive import ArchCompare

def main():
    
    usage = "\n %prog [options] -a archive_a.tar -b archive_b.tar -j fileTypes.json -c name"
    
    optParser = argparse.ArgumentParser()
    optParser.add_argument( "-a", "--archive_a", type=str, dest="archive_a",
    default = "", help = "archive path for a folder, file or a tar data" )
    
    optParser.add_argument( "-b", "--archive_b", type=str, dest="archive_b",
    default = "", help = "archive path for a folder, file or a tar data" )
    
    optParser.add_argument( "-j", "--json_config", type=str, dest="json_config",
    default = "", help = "path to json config file" )
    
    optParser.add_argument( "-c", "--cmp_type", type=str, dest="cmp_type", required=False, \
    default = None, help = "Compariosn type to perform [ \
     name: compares archives using file names, \
     checksum: perform comparsion based on ckecsum tool defined in json config file \
     data: does full data comparison based on tools defined for each extension type \
     Note- command line option if set overrides default value in json config file" )
    
    optParser.add_argument("-v", action="store_true", dest="verbose", default=False,help = "suppress progress report and warnings")
    optParser.add_argument("-q", action="store_false", dest="verbose", default=True)

    if len( sys.argv ) == 1:
        optParser.print_help()
        sys.exit( 1 )
    opts=optParser.parse_args()
    print('Running analysis, please check log file for any errors')
    mycomp=ArchCompare(archive_a=opts.archive_a,archive_b=opts.archive_b,json_config=opts.json_config,cmp_type=opts.cmp_type)
    mycomp.run_comparison()

if __name__ == '__main__':

  main()
