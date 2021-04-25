import argparse
import sys
from os import path
from pysodist.utils import utilities

log = utilities.log
vlog = utilities.vlog

''' This tool is designed to help configure pysodist for a given dataset. It is not necessary, but may make running
pysodist easier. To use it, you should first know the following about your dataset
    * MS 1 resolution - pysodist will use this value to make an initial guess at the Gaussian width used in fitting. 
    approximate peak width for your peptides (pysodist with extract this width +/- 
    * Approximate peak width - pysodist will use this as the expected value for the peak width
    * The residue labeling file you expect to use
    * The atom labeling file you expect to us
    * The path to the isodist executable you expect to use
'''


def add_args(parser):
    parser.add_argument('output', help='the configuration file you would like to output')
    parser.add_argument('--logfile', default=None, help='Optionally provide a path to a logfile to store outputs')
    return parser


def main(args):
    log('****INITIATING****', args.logfile)
    log('executed command: ' + " ".join(sys.argv), args.logfile)
    print('>> Please enter the instrument MS1 resolution at 200 m/z (e.g. QE-HFX: 60000; Sciex 5600: 25000.')
    ms1_res = float(input().replace(',', ''))
    print('>> Please enter the expected chromatographic peak width (FWHM) in seconds for a typical peptide (e.g. 10)')
    peak_width = float(input())
    print('>> Please enter the absolute path to the isodist labeling file you would like to use.')
    modelfile = input().replace('\\', '/')
    while not path.exists(modelfile):
        print('>> I could not find this file, please try entering the path again (and check for upper/lower case).')
        modelfile = input().replace('\\', '/')
    print('>> Please enter the absolute path to the isodist atom file you would like to use.')
    atomfile = input().replace('\\', '/')
    while not path.exists(atomfile):
        print('>> I could not find this file, please try entering the path again (and check for upper/lower case).')
        atomfile = input().replace('\\', '/')
    print('>> Please enter the absolute path to the isodist executable you would like to use.')
    isodist_executable = input().replace('\\', '/')
    while not path.exists(isodist_executable):
        print('>> I could not find this file, please try entering the path again (and check for upper/lower case).')
        isodist_executable = input().replace('\\', '/')
    with open(args.output, 'w') as f:
        f.write('FIELD\tVALUE\n')
        f.write('resolution\t' + str(ms1_res)+'\n')
        f.write('peak_width\t' + str(peak_width)+'\n')
        f.write('modelfile\t' + modelfile+'\n')
        f.write('atomfile\t' + atomfile+'\n')
        f.write('isodist_executable\t' + isodist_executable+'\n')
    log('\n++++COMPLETED configure++++', args.logfile)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description='Pysodist configuration. Used to generate a configuration file that pysodist can use for a given'
                    'dataset (or groups of related datasets). This step is not necessary, but may make running pysodist'
                    'easier. To use it, you should first know the following about your dataset: '
                    '   * MS 1 resolution - we use this in the initial guess of the Gaussian width used to fit peaks. '
                    '   * Approximate peak width - pysodist will use this as the expected value for the peak width. '
                    '   * The model labeling file you expect to use.'
                    '   * The atom labeling file you expect to use.'
                    '   * The path to the isodist executable you expect to use.')
    add_args(argparser)
    main(argparser.parse_args())
