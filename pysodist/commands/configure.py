# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
@version: 0.0.5
"""

import argparse
import sys
from os import path
from pysodist.utils import utilities
import pandas as pd
import pysodist

log = utilities.log

''' This tool is designed to help configure pysodist for a given dataset. To use it, you should must know the following
    * MS 1 resolution - pysodist will use this value to make an initial guess at the Gaussian width used in fitting. 
    approximate peak width for your peptides (pysodist with extract this width +/- 
    * Approximate peak width - pysodist will use this as the expected value for the peak width
    * The residue labeling file you expect to use
    * The atom labeling file you expect to us
    * The path to the isodist executable you expect to use
'''


def add_args(parser):
    parser.add_argument('output_directory', help='The base folder where all of our pysodist outputs will be saved.')
    parser.add_argument('--preconfigured', type=str, default=None,
                        help='Provide the full path to a configuration file. If provided, all other options will be'
                             'ignored, and this configuration file will simply be checked for errors/omissions.')
    parser.add_argument('--isodist_executable', type=str, default='PYTHON',
                        help='Full path to a compiled isodist executable if you want to use the Fortran implementation.'
                             ' Default is to use the python implementation, which is specified by "PYTHON" and does'
                             'not require providing a path.')
    parser.add_argument('--sample_list', nargs='*', default=None,
                        help='An optional list of samples to parse. By Default '
                             'all samples in the report are analyzed. Each sample separated by a space')
    parser.add_argument('--protein_list', nargs='*', default=None,
                        help='An optional list of the proteins to parse. By Default, all proteins in the report are '
                             'analyzed. Each Protein Gene Name separated by a space.')
    parser.add_argument('--isotope', type=str, default='light',
                        help='By Default, it is assumed that the report contains a light '
                             'isotope (no special labeling), if this field is not present'
                             'in the report, you can specify a different field here '
                             '(e.g. "heavy")')
    parser.add_argument('--q_value', type=float, default=0.00,
                        help='Used to optionally filter the report file based on the q_value. '
                             'By default, no q_value filtering is used.')
    # noinspection PyProtectedMember
    parser.add_argument('--atom_file', type=str, default=pysodist._ROOT+'/utils/model_files/atoms.txt',
                        help='Absolute path to the atom file '
                             '(typically in [pysodist_installed_directory]/model_files/atoms.txt')
    # noinspection PyProtectedMember
    parser.add_argument('--res_file', type=str, default=pysodist._ROOT+'/utils/model_files/U_var500N_fix998N.txt',
                        help='Absolute path to the residue file '
                             '(typically in [pysodist_installed_directory]/model_files/U_var500N_fix998N.txt)'
                             'Note that you may need to create a new file based on your labeling scheme. When you'
                             'actually run your fitting algorithm, you can also specify a new residue file if you want'
                             'to iteratively try different labeling schemes to best fit your data.')
    parser.add_argument('--ms1_resolution', type=float, default=60000,
                        help='MS1 resolution (calculated as full width half max at m/z=200).'
                             'Typical values are 25000 (SCIEX 5600); 60000 (QE HF-x). Default is 60000.')
    parser.add_argument('--peak_rt_width', type=float, default=10,
                        help='Typical peak width in seconds (calculated as full width half max).'
                             'Typical values are 5-20. Default is 10.')
    parser.add_argument('--logfile_name', default=None, help='Optionally provide a logfile name to store outputs')
    return parser


def main(args):
    print(args.atom_file)
    args.output_directory = args.output_directory.replace('\\', '/')
    if args.output_directory[-1] != '/':
        args.output_directory += '/'

    try:
        logfile = args.output_directory + args.logfile_name
    except TypeError:
        logfile = None
    log('****INITIATING****', logfile)
    log('executed command: ' + " ".join(sys.argv), logfile)
    if args.preconfigured is None:
        log('No pre-configuration file provided, checking arguments provided at the command line.')
        data = {'VALUE': {'isodist_executable': args.isodist_executable,
                          'sample_list': args.sample_list,
                          'protein_list': args.protein_list,
                          'isotope': args.isotope,
                          'atom_file': args.atom_file,
                          'res_file': args.res_file,
                          'q_value': args.q_value,
                          'ms1_resolution': args.ms1_resolution,
                          'peak_rt_width': args.peak_rt_width,
                          'logfile_name': args.logfile_name}}
        config_data = pd.DataFrame(data=data)
        config_data.index.name = 'FIELD'
    else:
        if path.exists(args.preconfigured):
            log('Pre-configuration file: ' + args.preconfigured + ' provided. Checking each argument.', logfile)
            config_data = pd.read_csv(args.preconfigured, sep=',', index_col='FIELD', comment='#')
        else:
            log('***ERROR*** Pre-configuration file: ' + args.preconfigured +
                ' was not found. Please check that the file exists and try again.', logfile)
            sys.exit('***ERROR*** Please check the log.')
    config_data.loc['output_directory'] = args.output_directory
    config_data.loc['isodist_executable']['VALUE'] = config_data.loc['isodist_executable']['VALUE'].replace('\\', '/')
    config_data.loc['atom_file']['VALUE'] = config_data.loc['atom_file']['VALUE'].replace('\\', '/')
    config_data.loc['res_file']['VALUE'] = config_data.loc['res_file']['VALUE'].replace('\\', '/')

    assert (config_data.loc['isodist_executable']['VALUE'] == 'PYTHON' or
            path.exists(config_data.loc['isodist_executable']['VALUE'])), \
        config_data.loc['isodist_executable']['VALUE'] + \
        ' not found. If planning to use the Python implementation, this field should be "PYTHON"'
    assert path.exists(config_data.loc['atom_file']['VALUE']), \
        log(config_data.loc['atom_file']['VALUE'] + ' not found.', logfile)
    assert path.exists(config_data.loc['res_file']['VALUE']), \
        log(config_data.loc['res_file']['VALUE'] + ' not found.', logfile)
    assert (0.0 <= float(config_data.loc['q_value']['VALUE']) <= 1.0)
    assert (1000 < float(config_data.loc['ms1_resolution']['VALUE']) < 1000000)
    assert (0 < float(config_data.loc['peak_rt_width']['VALUE']) < 200)

    log('all inputs were valid.', logfile)
    log('writing output configuration file to the output directory: ' +
        config_data.loc['output_directory']['VALUE'] + '00_config.cfg', logfile)
    with open(config_data.loc['output_directory']['VALUE'] + '00_config.cfg', 'w') as f:
        f.write('#Pysodist configuration file v0.0.5\n'
                '#Each line corresponds to a parameter, with the parameter name first, and the value second.\n'
                '#Field names and values should be separated by a comma.\n'
                '#The first used line should contain "FIELD,VALUE"\n')
    config_data.to_csv(config_data.loc['output_directory']['VALUE'] + '00_config.cfg', mode='a')
    log('\n++++COMPLETED configure++++', logfile)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description='Pysodist configuration. Used to generate a configuration file that pysodist can use for a given'
                    'dataset (or groups of related datasets). '
                    'To use it, you should first know the following about your dataset: '
                    '   * MS 1 resolution - we use this in the initial guess of the Gaussian width used to fit peaks. '
                    '   * Approximate peak width - pysodist will use this as the expected value for the peak width. '
                    '   * The model labeling file you expect to use.'
                    '   * The atom labeling file you expect to use.'
                    '   * The path to the isodist executable you expect to use.')
    add_args(argparser)
    main(argparser.parse_args())
