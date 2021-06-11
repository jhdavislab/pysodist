# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
@version: 0.0.5
"""

import argparse
import sys
from os import path, mkdir
from pysodist.utils import utilities
import pandas as pd
import pysodist

log = utilities.log
clean_path = utilities.clean_path

''' This tool is designed to help configure pysodist for a given dataset. To use it, you should must know the following
    * MS 1 resolution - pysodist will use this value to make an initial guess at the Gaussian width used in fitting. 
    approximate peak width for your peptides (pysodist with extract this width +/- 
    * Approximate peak width - pysodist will use this as the expected value for the peak width
    * The residue labeling file you expect to use
    * The atom labeling file you expect to us
    * The path to the isodist executable you expect to use
'''


def load_config_file(config_file, logfile_name=None):
    if not path.exists(config_file):
        sys.exit('***ERROR*** Pre-configuration file: ' + config_file +
                 ' was not found. Please check that the file exists and try again.')
    else:
        config_data = pd.read_csv(config_file, sep=',', index_col='FIELD', comment='#')
        config_data.loc['output_directory']['VALUE'] = clean_path(config_data.loc['output_directory']['VALUE'])

        if logfile_name is None:
            try:
                logfile = config_data.loc['output_directory']['VALUE'] + config_data['logfile_name']
            except TypeError:
                logfile = None
        else:
            logfile = config_data.loc['output_directory']['VALUE'] + logfile_name

        log('Pre-configuration file: ' + config_file + ' provided. Checking each argument.', logfile)

        config_data.loc['isodist_exe']['VALUE'] = config_data.loc['isodist_exe']['VALUE'].replace('\\', '/')
        config_data.loc['atom_file']['VALUE'] = config_data.loc['atom_file']['VALUE'].replace('\\', '/')
        config_data.loc['res_file']['VALUE'] = config_data.loc['res_file']['VALUE'].replace('\\', '/')

        assert (config_data.loc['isodist_exe']['VALUE'] == 'PYTHON' or
                path.exists(config_data.loc['isodist_exe']['VALUE'])), \
            config_data.loc['isodist_exe']['VALUE'] + \
            ' not found. If planning to use the Python implementation, this field should be "PYTHON"'
        assert path.exists(config_data.loc['atom_file']['VALUE']), \
            log(config_data.loc['atom_file']['VALUE'] + ' not found.', logfile)
        assert path.exists(config_data.loc['res_file']['VALUE']), \
            log(config_data.loc['res_file']['VALUE'] + ' not found.', logfile)
        assert (0.0 <= float(config_data.loc['q_value']['VALUE']) <= 1.0)
        assert (1000 < float(config_data.loc['ms1_resolution']['VALUE']) < 1000000)
        assert (0 < float(config_data.loc['peak_rt_width']['VALUE']) < 200)
        log('configure stage inputs were valid.', logfile)
    return config_data


def write_config_file(config_data, config_file_name, logfile=None):
    config_file = config_data.loc['output_directory']['VALUE'] + config_file_name
    log('writing output configuration file to the output directory: ' + config_file, logfile)
    with open(config_file, 'w') as f:
        f.write('# Pysodist configuration file v0.0.5\n'
                '# Each line corresponds to a parameter, with the parameter name first, and the value second.\n'
                '# Field names and values should be separated by a comma.\n'
                '# The first used line should contain "FIELD,VALUE"\n')
    config_data.to_csv(config_file, mode='a')


def add_args(parser):
    parser.add_argument('output_directory', help='The base folder where all of our pysodist outputs will be saved.'
                                                 'If a preconfigured file is provided, that output directory will'
                                                 'be used instead of the one provided here.')
    parser.add_argument('--preconfigured', type=str, default=None,
                        help='Provide the full path to a configuration file. If provided, all other options will be'
                             'ignored, and this configuration file will simply be checked for errors/omissions.')
    parser.add_argument('--logfile_name', default=None, help='Optionally provide a logfile name to store outputs')

    pi_group = parser.add_argument_group('parse_input')
    es_group = parser.add_argument_group('extract_spectra')
    fs_group = parser.add_argument_group('fit_spectra')
    af_group = parser.add_argument_group('analyze_fit')
    pi_group.add_argument('--guide_file', type=str, default=None,
                          help='Provide the full path to the guide file (e.g. a skyline report file). If not provided, '
                               'the default of None is used and this file can be provided during the parse_input.')
    pi_group.add_argument('--sample_list', nargs='*', default=None,
                          help='An optional list of samples to parse. By Default '
                               'all samples in the report are analyzed. Each sample separated by a space')
    pi_group.add_argument('--protein_list', nargs='*', default=None,
                          help='An optional list of the proteins to parse. By Default, all proteins in the report are '
                               'analyzed. Each Protein Gene Name separated by a space.')
    pi_group.add_argument('--isotope', type=str, default='light',
                          help='By Default, it is assumed that the report contains a light '
                               'isotope (no special labeling), if this field is not present'
                               'in the report, you can specify a different field here '
                               '(e.g. "heavy")')
    pi_group.add_argument('--q_value', type=float, default=0.00,
                          help='Used to optionally filter the report file based on the q_value. '
                               'By default, no q_value filtering is used.')

    es_group.add_argument('--mzml_file', type=str, default=None,
                          help='Provide the full path to the mzml file. If not provided, the default of None is used, '
                               'and this file can be provided during the extract_spectra stage.')

    fs_group.add_argument('--isodist_exe', type=str, default='PYTHON',
                          help='Full path to a compiled isodist executable if you want to use the Fortran version.'
                               ' Default is to use the python implementation, which is specified by "PYTHON" and does'
                               'not require providing a path.')
    # noinspection PyProtectedMember
    fs_group.add_argument('--atom_file', type=str, default=pysodist._ROOT + '/model_files/atoms.txt',
                          help='Absolute path to the atom file '
                               '(typically in [pysodist_installed_directory]/model_files/atoms.txt')
    # noinspection PyProtectedMember
    fs_group.add_argument('--res_file', type=str, default=pysodist._ROOT + '/model_files/U_var500N_fix998N.txt',
                          help='Absolute path to the residue file '
                               '(typically in [pysodist_installed_directory]/model_files/U_var500N_fix998N.txt)'
                               'Note that you may need to create a new file based on your labeling scheme. When you'
                               'actually run fit_spectra, you can also specify a new residue file if you want'
                               'to iteratively try different labeling schemes to best fit your data.')
    fs_group.add_argument('--ms1_resolution', type=float, default=60000.0,
                          help='MS1 resolution (calculated as full width half max at m/z=200).'
                               'Typical values are 25000 (SCIEX 5600); 60000 (QE HF-x). Default is 60000.')

    af_group.add_argument('--peak_rt_width', type=float, default=10.0,
                          help='Typical peak width in seconds (calculated as full width half max).'
                               'Typical values are 5-20. Default is 10.')
    return parser


def main(args):
    output_directory = clean_path(args.output_directory)
    if path.exists(output_directory):
        print('output directory already exists, using this directory.')
    else:
        print('output directory does not exist. Creating directory: ' + output_directory)
        mkdir(output_directory)

    try:
        logfile = output_directory + args.logfile_name
    except TypeError:
        logfile = None

    log('++++INITIATING++++', logfile)
    log('executed command: ' + " ".join(sys.argv), logfile)
    if args.preconfigured is None:
        log('no pre-configuration file provided, checking arguments provided at the command line.', logfile)
        data = {'VALUE': {'isodist_exe': args.isodist_exe,
                          'sample_list': args.sample_list,
                          'protein_list': args.protein_list,
                          'isotope': args.isotope,
                          'atom_file': args.atom_file,
                          'res_file': args.res_file,
                          'q_value': args.q_value,
                          'ms1_resolution': args.ms1_resolution,
                          'peak_rt_width': args.peak_rt_width,
                          'logfile_name': args.logfile_name,
                          'output_directory': args.output_directory,
                          'guide_file': args.guide_file,
                          'mzml_file': args.mzml_file}}
        config_data = pd.DataFrame(data=data)
        config_data.index.name = 'FIELD'
        config_data = config_data.reindex(['output_directory',
                                           'guide_file',
                                           'mzml_file',
                                           'sample_list',
                                           'protein_list',
                                           'q_value',
                                           'isotope',
                                           'isodist_exe',
                                           'atom_file',
                                           'res_file',
                                           'ms1_resolution',
                                           'peak_rt_width',
                                           'logfile_name'])
    else:
        config_data = load_config_file(args.preconfigured, logfile)

    write_config_file(config_data, '00_config.cfg', logfile)
    log('++++COMPLETED configure++++', logfile)


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
