# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
@version: 0.0.4
"""

import pysodist.commands.parse_input as parse_input
import pysodist.commands.extract_spectra as extract_spectra
import pysodist.commands.run_isodist as run_isodist
import pysodist.commands.plot_spectra as plot_spectra
import pysodist
import os
import argparse
import pandas as pd
import math
import subprocess


def add_args(parser):
    # required
    parser.add_argument('input', help='input file to parse. Currently only skyline report files are supported')
    parser.add_argument('mzml', help='the relative path to mzml file to be analyzed. For Thermo instruments, '
                                     'one should generate the .mzml file from the original .raw file using msconvert '
                                     'as follows: \ .\msconvert.exe ".\[YOUR_RAW_FILE].raw" -o "./" --mzML --64 -v '
                                     'mz64 --inten32 --noindex --filter "msLevel 1" --zlib')
    parser.add_argument('sample_name', help='name of the sample within the skyline report to be analyzed.')
    parser.add_argument('isodist_command', help='exact fortran command to execute. e.g. C:\isodist\isodist.exe')
    parser.add_argument('atomfile',
                        help='Specify the path to the atom definition file (e.g. exp_atom_defs.txt). You will likely '
                             'not need to modify this file.')
    parser.add_argument('resfile', help='Specify the path to the residue labeling file - you will likely need to edit '
                                        'this file based on your labeling scheme.\ Note that your output will use the '
                                        'name of this file to ensure you know which model file produced which output.')

    # fortran isodist specific parser.add_argument('--batch_size', default=50, type=int, help='Number of peptides to
    # be fit in each batch by isodist.')
    parser.add_argument('--threads', default=4, type=int,
                        help='number of threads to use. typically 1 less than the number of cores available. Default=4')
    parser.add_argument('--wait_time', default=60, type=int,
                        help='number of seconds to wait between each polling to test if the isodist run has finished. '
                             'Default=60 seconds')

    # extract spectra specific
    parser.add_argument('--output_directory', default='./',
                        help='Output files will be saved in this folder: 1 directory per sample in the skyline '
                             'report. Default = ./')
    parser.add_argument('--protein_list', default=None,
                        help='An optional list of the proteins to parse. By default, all proteins in the report are '
                             'analyzed.')
    parser.add_argument('--isotope', default='light',
                        help='Be default, it is assumed that the report contains a light isotope (no special '
                             'labeling), if this field is not present in the report, you can specify a different '
                             'field here (e.g. "heavy")')
    parser.add_argument('--q_value', default=0.00,
                        help='Used to optionally filter the report file based on the q_value. By default, no q_value '
                             'filtering is used.')
    parser.add_argument('--labeling', default='N15',
                        help='The labeling scheme used for the highest mass isotope envelope you expect to fit. E.g. '
                             'N15 or C13')
    parser.add_argument('--interp_only', action='store_const', const=True, default=False,
                        help='Only save the interpolated spectra instead of the raw spectra')
    parser.add_argument('--sum_only', action='store_const', const=True, default=False,
                        help='Only save summed (and interpolated) spectra instead of all individual spectra. Results '
                             'in 1 spectra per peptide')
    parser.add_argument('--interp_res', default=0.001, type=float,
                        help='Set the interpolation delta m/z - typical values from 0.01 to 0.001')

    # plotting specific
    parser.add_argument('--numerator', nargs='+', default=['AMP_U'],
                        help='list of the fields to use in the numerator of the abundance ratio calculation ('
                             'typically AMP_U, AMP_L, AMP_F, or some combination. Default is AMP_U.')
    parser.add_argument('--denominator', nargs='+', default=['AMP_U', 'AMP_F'],
                        help='list of the fields to use in the denominator of the abundance ratio calculation ('
                             'typically AMP_U, AMP_L, AMP_F, or some combination. Default is AMP_U, AMP_F')
    parser.add_argument('--no_png', action='store_const', const=True, default=False,
                        help='By default .png files for the plots will be saved. This option forces these to not be '
                             'saved.')
    parser.add_argument('--no_pdf', action='store_const', const=True, default=False,
                        help='By default .pdf files for the plots will be saved. This option forces these to not be '
                             'saved.')

    return parser


def main(args):
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('***stage 1: parsing input file...***')
    sample_name = args.sample_name.strip()
    isodist_df = parse_input.parse_skyline(args.input, protein_list=args.protein_list, sample_list=[sample_name],
                                           isotope=args.isotope, q_value=args.q_value,
                                           output_directory=args.output_directory, IO=True)[0]
    print('unique peptides found: ' + str(isodist_df.shape[0]))
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('***stage 2: extracting spectra...***')
    output_directory = args.output_directory.replace('\\', '/')
    if output_directory[-1] != '/':
        output_directory += '/'
    parsed_report = output_directory + args.sample_name + '/pd_parsed_report.tsv'
    print('using pysodist report file: ' + parsed_report)
    assert (os.path.exists(parsed_report) is True)
    parsed_mzml = extract_spectra.parse_mzml(args.mzml)
    sample_output_directory = output_directory + args.sample_name + '/'
    assert (os.path.exists(sample_output_directory) is True)
    extract_spectra.extract_spectra(parsed_mzml, parsed_report, sample_output_directory,
                                    labeling=args.labeling, save_interp_spectra=args.interp_only,
                                    interp_res=args.interp_res, sum_spectra_only=args.sum_only)
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('***stage 3: fitting spectra using isodist...***')
    isodist_input_file = sample_output_directory + 'pd_exported_peaks.tsv'
    atomfile = args.atomfile.replace('\\', '/')
    resfile = args.resfile.replace('\\', '/')
    isodist_command = args.isodist_command.replace('\\', '/')
    resfile_name = resfile.split('/')[-1].split('.txt')[0]
    isodist_output_csv = sample_output_directory + resfile_name + '_output.csv'

    print('working in directory: ' + sample_output_directory)
    run_isodist.prep_model_files(sample_output_directory, atomfile, resfile)

    batch_base_path = '/'.join(isodist_input_file.split('/')[:-1]) + '/'
    batch_df = pd.read_csv(isodist_input_file, sep='\t')
    num_spectra = batch_df.shape[0]
    batch_size = math.ceil(num_spectra / args.threads)
    batch_file_path_list = run_isodist.write_batch_files(batch_df, batch_base_path, batch_size=batch_size)

    in_file_list = []
    for batch_file_path in batch_file_path_list:
        batch_file_path = batch_file_path.replace('\\', '/')
        in_file_list.append(run_isodist.write_isodist_input(batch_file_path, atomfile, resfile))

    csv_list = run_isodist.run_fortran_isodist(in_file_list, isodist_command, threads=args.threads,
                                               wait_time=args.wait_time)
    run_isodist.compile_isodist_csvs(csv_list, isodist_output_csv, parsed_pysodist_input=parsed_report)
    print('cleaning up...')
    run_isodist.cleanup(isodist_output_csv)
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('***stage 4: generating plots of the results...***')

    isodist_compiled_csv_name = isodist_output_csv.split('/')[-1]
    isodist_csv_no_extension = isodist_compiled_csv_name.split('.csv')[0]
    base_folder_name = isodist_output_csv.split('_output.csv')[0]
    fit_folder = base_folder_name + '_isodist_fits/'
    isodist_output_folder = base_folder_name + '_isodist_outputs/'
    isodist_output_csv = isodist_output_folder + isodist_compiled_csv_name
    print('parsing isodist csv file: ' + isodist_output_csv)
    isodist_output_pd = plot_spectra.parse_isodist_csv(isodist_output_csv)

    current_ratio_string = plot_spectra.get_current_ratio_string(args.numerator, args.denominator, isodist_output_pd)
    print('all of the following plots will use current ratio as: ' + current_ratio_string)
    isodist_output_pd = plot_spectra.set_current_ratio(isodist_output_pd, numerator=args.numerator,
                                                       denominator=args.denominator)

    plot_output_folder = base_folder_name + '_final_plots/'
    try:
        os.mkdir(plot_output_folder)
    except OSError:
        print(
            '...the output directory: ' + plot_output_folder + 'already exists, and files within it may be '
                                                               'overwritten. continue? [y/n]')
        choice = input().lower()
        if not choice == 'y':
            raise

    isodist_output_pd.to_csv(plot_output_folder + isodist_csv_no_extension + '_isodist_result.csv')

    assert not (args.no_png and args.no_pdf)

    print('using fit spectra from directory: ' + fit_folder)
    all_proteins = list(set(isodist_output_pd['protein'].values))

    print('saving plots for ' + str(len(all_proteins)) + ' proteins to: ' + plot_output_folder)
    all_proteins.sort()
    for protein in all_proteins:
        print('plotting spectra for protein: ' + protein)
        related_spectra = plot_spectra.get_by_protein(isodist_output_pd, protein)
        plot_spectra.plot_spectra_group(related_spectra, fit_folder, working_path=sample_output_directory,
                                        numerator=args.numerator,
                                        denominator=args.denominator, png=not args.no_png,
                                        pdf=not args.no_pdf, saved_output=plot_output_folder + protein)

    print('calculating and plotting abundances...')
    plot_spectra.plot_all_ratios(isodist_output_pd, numerator=args.numerator, denominator=args.denominator,
                                 saved_output_path=plot_output_folder, png=not args.no_png, pdf=not args.no_pdf)
    print('plotting the csv stat histograms...')
    plot_spectra.plot_csv_stats(isodist_output_pd, current_ratio_string, output_path=plot_output_folder,
                                png=not args.no_png, pdf=not args.no_pdf)

    print('copying analysis_template.ipynb jupyter notebook...')
    out_ipynb = plot_output_folder + 'pysodist_analysis.ipynb'
    if not os.path.exists(out_ipynb):
        # noinspection PyProtectedMember
        root_path = pysodist._ROOT + '/'
        ipynb = root_path + 'utils/analysis_template.ipynb'
        cmd = f'cp {ipynb} {out_ipynb}'
        subprocess.check_call(cmd, shell=True)
    else:
        print(f'{out_ipynb} already exists. Skipping')


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description='Pysodist runner - links together the various pysodist modules (parse_input, extract_spectra, '
                    'run_isodist, plot_spectra. Note that this is experimental and should only be used if you carefully'
                    'inspect the code for errors.')
    add_args(argparser)
    main(argparser.parse_args())
