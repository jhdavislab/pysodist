# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
"""

import numpy as np
import pandas as pd
import pickle
import pysodist.utils.skyline_report_defs as defs
from pyteomics import mzml
import re
import os
import argparse
from scipy.interpolate import interp1d
from pysodist.utils import utilities
import sys

log = utilities.log
vlog = utilities.vlog


def write_scan(select_scan, file_name):
    """
    Writes an individual scan to a file that can be read by default (fortran) isodist
    :param select_scan: a pandas dataframe with the selected scan (only including the desired m/z range).
    Can be interpolated, or summed, or raw.
    :param file_name: string pointing to the full path for the file to be written.

    returns: None
    """
    select_scan.to_csv(file_name, sep=' ', header=False, index=False)
    return None


def extract_spectra(parsed_mzml, parsed_report, output_dir,
                    labeling='N15', save_interp_spectra=False, interp_res=0.001, sum_spectra_only=False, logfile=None):
    """
    Extracts spectra as specified in the parsed level 1 output. Writes individual spectra files as well as
    a batch file for them. Optionally can interpolate the spectra, and sum them on a per-peptide basis.

    :param parsed_mzml: the parsed mzml file (should be output of parse_mzml)
    :param parsed_report: the full path to output file from the level1 parser.
    :param output_dir: the full path to the 'base' directory. Within it, a spectra will be created.
    :param labeling: string defining the labeling type, which is used to determine the width of the m/z window.
                     Currently implemented N15 and C13.
    :param save_interp_spectra: bool defining whether to save interpolated spectra instead of raw spectra (spectra are
    interpolated at interp_res)
    :param interp_res: the resolution to interpolate to when producing the summed spectra. Default 0.001.
    :param sum_spectra_only: boolean determining if only the summed spectra should be saved
    (single spectra per peptide). Uses the interp_res to do the summation. Default False.
    :param logfile: path to file to log

    returns: pandas dataframe ready to be written as a batch file(s).
    """

    parsed_report_df = pd.read_csv(parsed_report, sep='\t')
    rt_array = parsed_mzml['retention_times']
    log('total spectra in mzml: ' + str(parsed_mzml['retention_times'].shape[0]), logfile)
    log('rt_range in mzml: ' + str(parsed_mzml['retention_times'][0]) + ' - ' +
        str(parsed_mzml['retention_times'][-1]) + ' mins', logfile)
    log('total peptides to extract: ' + str(parsed_report_df.shape[0]), logfile)
    log('peptide rt_range: ' + str(parsed_report_df['rt_start'].min()) + ' - ' +
        str(parsed_report_df['rt_end'].max()) + ' mins', logfile)
    spectra_dict = {}

    peaks_dir = output_dir + 'spectra/'
    local_peaks_dir = './spectra/'
    try:
        os.mkdir(peaks_dir)
    except OSError:
        print('...the output spectra directory: ' + peaks_dir +
              'already exists, and files within it may be overwritten. continue? [y/n]')
        choice = input().lower()
        if not choice == 'y':
            raise

    for index, current_peptide in parsed_report_df.iterrows():
        first_scan = np.argmax(rt_array >= current_peptide['rt_start'])
        last_scan = np.argmax(rt_array > current_peptide['rt_end']) - 1
        mz_range = mz_window(current_peptide['peptide_modified_sequence'],
                             current_peptide['charge'], current_peptide['mz'], labeling=labeling)

        if (last_scan - first_scan) < 2:
            log('peptide ' + str(index) + ' : ' + current_peptide['peptide_modified_sequence'] +
                'has a very narrow RT range with only ' + str(last_scan-first_scan+1) +
                ' scans. It is being skipped.', logfile)
        else:
            log('extracting ' + str(last_scan - first_scan + 1) + ' spectra for peptide ' + str(index) + ' : ' +
                current_peptide['peptide_modified_sequence'], logfile)

            interp_mz_axis = np.arange(mz_range[0], mz_range[1], interp_res)
            interp_summed_intensity = np.zeros(interp_mz_axis.shape[0])

            peptide_base_name = "_".join([current_peptide['peptide_modified_sequence'],
                                          str(current_peptide['charge']),
                                          str(round(current_peptide['mz'], 3)),
                                          str(current_peptide['start_pos']) + '-' + str(current_peptide['end_pos'])])

            for current_scan_num in range(first_scan, last_scan+1):
                scan_data = parsed_mzml['ms1_scans'][current_scan_num]
                select_scan_data = scan_data.loc[
                    (scan_data['mz_data'] >= mz_range[0]) & (scan_data['mz_data'] <= mz_range[1])]
                if len(select_scan_data) < 4:
                    log('**++** PEPTIDE: ' + str(current_peptide['peptide_modified_sequence']) + ', scan #: ' + str(
                        current_scan_num) + 'has fewer than 4 MS1 points in  m/z range. Skipping this scan.', logfile)
                else:
                    interpolation = interp1d(scan_data['mz_data'], scan_data['intensity_data'])
                    interp_intensity = interpolation(interp_mz_axis)
                    interp_scan_data = pd.DataFrame({'mz_data': interp_mz_axis, 'intensity_data': interp_intensity})

                    interp_summed_intensity += interp_intensity

                    if sum_spectra_only is False:
                        file_name = peptide_base_name + '_' + str(round(rt_array[current_scan_num], 3))
                        spectra_string = peaks_dir + file_name + '.tsv'
                        local_spectra_string = local_peaks_dir + file_name + '.tsv'
                        if save_interp_spectra:
                            write_scan(interp_scan_data, spectra_string)
                        else:
                            write_scan(select_scan_data, spectra_string)

                        spectra_dict[file_name] = [current_peptide['peptide_modified_sequence'],
                                                   str(current_peptide['charge']),
                                                   local_spectra_string]

            file_name = peptide_base_name + '_SUM'
            spectra_string = peaks_dir + file_name + '.tsv'
            local_spectra_string = local_peaks_dir + file_name + '.tsv'
            interp_summed_spectra = pd.DataFrame(
                {'mz_data': interp_mz_axis, 'intensity_data': interp_summed_intensity / ((last_scan - first_scan) / 2)})
            write_scan(interp_summed_spectra, spectra_string)
            spectra_dict[file_name] = [current_peptide['peptide_modified_sequence'],
                                       str(current_peptide['charge']),
                                       local_spectra_string]

            spectra_df = pd.DataFrame.from_dict(spectra_dict, orient='index',
                                                columns=['peptide_modified_sequence', 'charge', 'spectra_file'])
            spectra_df.to_csv(output_dir + 'pd_exported_peaks.tsv', sep='\t', index=False)
    return pd.DataFrame.from_dict(spectra_dict, orient='index',
                                  columns=['peptide_modified_sequence', 'charge', 'spectra_file'])


def parse_mzml(mzml_path, pickle_data=None, logfile=None):
    """
    retrieves all scans from a portion of an mzml file and generates arrays of mz versus intensity at each RT

    :param mzml_path: string pointing to the .mzml file to extra scans from
    :param pickle_data: string pointing to a pickle file to save the extracted spectra into. Default is not saved
    :param logfile: path to file to log

    :return: a dictionary with keys retention_times, and ms1_scans. retention_times points to a numpy array of
    each retention time. ms1_scans points to a list of pandas dataframes with columns mz_data and intensity_data.
    The index of this list corresponds to the position in the retention_times numpy array.
    """
    with mzml.read(mzml_path) as mz_reader:
        scan_list = []
        rt_list = []
        log('reading mzml file ' + mzml_path + '...', logfile)
        for scan in mz_reader:
            assert scan['ms level'] == 1, 'Your mzml file contains non-MS1 scans. ' \
                                          'When converting your mzml file, only include MS1 scans. See the pysodist' \
                                          'docs for how to appropriately convert to .mzml files using msconvert.'
            mz_int_pd = pd.DataFrame({'mz_data': scan['m/z array'], 'intensity_data': scan['intensity array']})
            scan_list.append(mz_int_pd)
            new_rt = float(scan['scanList']['scan'][0]['scan start time'])
            rt_list.append(new_rt)
    parsed_mz_file = {'retention_times': np.array(rt_list), 'ms1_scans': scan_list}
    if not (pickle_data is None):
        pickle.dump(parsed_mz_file, open(pickle_data, 'wb'))
    return parsed_mz_file


def mz_window(peptide_modified_sequence, z, mz, labeling='N15', topomers_left=2, tomopers_right=5):
    """
    calculates the mz window required to accommodate the range of mzs  for all isotopes present in a given peptide

    :param peptide_modified_sequence: string of peptide amino acid sequence including modifications
    :param z: int of peptide charge state
    :param mz: float of mz for the unlabeled monoisotopic peptide
    :param labeling: isotope labeling method ie. "K8R10" or "N15" (default None: unlabeled)
    :param topomers_left: number of isotopomers to include in the extracted spectra to the left of the min mass
    :param tomopers_right: number of isotopomers to include in the extracted spectra to the right of the max mass

    :return: list of lower and upper range of possible mz values +/- 8 mass units for each peptide
    """
    peptide_simple_seq = ''.join(re.split('\[\+\d+\.\d+]', peptide_modified_sequence))
    peptide_simple_seq = peptide_simple_seq.upper()
    n14_offset = defs.N15MASS - defs.N14MASS
    c13_offset = defs.C13MASS - defs.C12MASS
    lys_count = peptide_simple_seq.count('K')
    arg_count = peptide_simple_seq.count('R')

    if labeling == 'N15':
        labeled_mz = ((sum([defs.AANITROGENS[i] * n14_offset for i in peptide_simple_seq])) / z) + mz
    elif labeling == 'C13':
        labeled_mz = ((sum([defs.AACARBONS[i] * c13_offset for i in peptide_simple_seq])) / z) + mz
    elif labeling == 'K8R10':
        labeled_mz = ((lys_count * 8 + arg_count * 10) / z) + mz
    elif labeling == 'K6R6':
        labeled_mz = ((lys_count * 6 + arg_count * 6) / z) + mz
    else:
        labeled_mz = mz
    return [mz - (topomers_left / z), labeled_mz + (tomopers_right / z)]


def add_args(parser):
    parser.add_argument('mzml', help='the relative path to mzml file to be analyzed.')
    parser.add_argument('parsed_report', help='the parsed report file generated by pysodist parse_input.py')
    parser.add_argument('--labeling', default='N15', help='The labeling scheme used for the highest '
                                                          'mass isotope envelope you expect to fit. '
                                                          'Possible values are: N15, C13, K6R6, K8R10.')
    parser.add_argument('--interp_only', action='store_const', const=True, default=False,
                        help='Only save the interpolated spectra instead of the raw spectra. Optional, default is N15.')
    parser.add_argument('--sum_only', action='store_const', const=True, default=False,
                        help='Only save summed (and interpolated) spectra instead of all individual spectra. '
                             'Results in 1 spectra per peptide. Optional, default is to extract all spectra.')
    parser.add_argument('--interp_res', default=0.001, type=float,
                        help='Set the interpolation delta m/z - typical values from 0.01 to 0.001. '
                             'Optional, default is 0.001.')
    parser.add_argument('--logfile', default=None, help='Optionally provide a path to a logfile to store outputs')
    return parser


def main(args):
    log('\n****INITIATING****', args.logfile)
    log('executed command: ' + " ".join(sys.argv), args.logfile)
    parsed_report = args.parsed_report.replace('\\', '/')
    output_dir = '/'.join(parsed_report.split('/')[:-1]) + '/'
    assert (os.path.exists(parsed_report) is True)

    parsed_mzml = parse_mzml(args.mzml, logfile=args.logfile)
    extract_spectra(parsed_mzml, args.parsed_report, output_dir,
                    labeling=args.labeling, save_interp_spectra=args.interp_only,
                    interp_res=args.interp_res, sum_spectra_only=args.sum_only)
    log('\n++++COMPLETED extract_spectra++++\n\n', args.logfile)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description='Pysodist .mzml parser, used to extract the relevant spectra for a single sample as defined by the'
                    'parsed report file (e.g. from Skyline or Encyclopedia) from the matching .mzml file. For Thermo'
                    'instruments, one should generate the .mzml file from the original .raw file using msconvert'
                    'as follows:.\msconvert.exe ".\[YOUR_RAW_FILE].raw" -o "./" --mzML --64 -v mz64 '
                    '--inten32 --noindex --filter "msLevel 1" --zlib. '
                    'This tool will extra the individual spectra as .tsv files (saved in a ./OUTPUT_peaks directory)'
                    'and product .batch and .in files to be used by Fortran isodist')
    add_args(argparser)
    main(argparser.parse_args())
