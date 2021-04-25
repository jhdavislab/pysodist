# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
"""

import pandas as pd
import pysodist.utils.skyline_report_defs as defs
import argparse
import os
import sys
from pysodist.utils import utilities

log = utilities.log
vlog = utilities.vlog

''' This tool is designed as an initial parser for skyline report files to get them into a common form that pysodist 
can then use to extract spectra. This has been designed to allow for rapid changes to the skyline parser if the report 
files change, and to allow other types of 'guides' (anything that tells you retention times to extract for given 
spectra, e.g. a MS/MS search) to be read and converted to a fixed file format the pysodist can handle. This common file 
format requires the following fields for each peptide:
    * rt_start (starting RT to extract MS1 from)
    * rt_end (ending RT to extract MS1 from)
    * peptide_modified_sequence (this should include the full monoisotopic mod names e.g. C[+57.021464])
    * charge
    * mz (this is for the unlabeled species)
    * protein_IDs (for skyline, this will typically be the Protein Gene Name)
    * peptide_start_position (AA number of the peptide start)
    * peptide_end_position (AA number of the peptide end)
'''


def map_int(i):
    """
    helper function to return ints instead of strings/floats
    """

    try:
        x = int(i)
    except ValueError:
        x = 0
    return x


def extract_skyline_sub(full_dataframe, sample, protein_list=None, isotope='light'):
    """
    extract_skyline_sub takes a full skyline report as input and outputs a filtered dataframe with only those columns
    relevant to the proteins and sample of interest

    :param full_dataframe: pandas dataframe directly uploaded from skyline
    :param sample: string of injection sample name
    :param protein_list: a list of strings of uniprot identifier for particular proteins of interest
    (default None includes all proteins)
    :param isotope: string for filtering output dataframe such that columns from a particular isotope type are included
    (default 'light')

    :return: a pandas dataframe with the only the sample and proteins of interest and only those columns needed for
    subsequent analyses
    """
    required_cols = [isotope + ' ' + sample + ' ' + defs.RT_START_FIELD,
                     isotope + ' ' + sample + ' ' + defs.RT_END_FIELD,
                     defs.PROTEIN_GENE_FIELD,
                     defs.PEPTIDE_MOD_SEQ_FIELD,
                     isotope + ' ' + defs.MZ_FIELD,
                     defs.PEPTIDE_CHARGE_FIELD,
                     isotope + ' ' + sample + ' ' + defs.DETECT_Q_VALUE_FIELD,
                     defs.START_POS_FIELD,
                     defs.END_POS_FIELD]  # fill in with the columns you need - use the isotope info to get correct ones
    new_data_frame = full_dataframe.reindex(columns=required_cols)
    assert required_cols[
               0] in full_dataframe.columns, 'the provided sample is likely not present. Check the sample name.'
    if not (protein_list is None):
        new_data_frame = new_data_frame[new_data_frame[defs.PROTEIN_GENE_FIELD].isin(protein_list)]
    return new_data_frame


def rename_peptide(skyline_peptide_modified_seq):
    isodist_peptide_modified_seq = skyline_peptide_modified_seq
    for pms_res, id_res in defs.ISODIST_MOD.items():
        isodist_peptide_modified_seq = isodist_peptide_modified_seq.replace(pms_res, id_res)
    return isodist_peptide_modified_seq


def parse_sub_skyline(skyline_sub_df, sample, q_value=0.00, isotope='light', logfile=None):
    """
    parse_sub_skyline is a helper function to parse a sample and protein filtered skyline dataframe such that a properly
    sorted and labeled dataframe is output

    :param skyline_sub_df: pandas dataframe with only one sample injection and only those proteins of interest
    :param sample: string of injection sample name
    :param q_value: float with a optional q_value cutoff to filter your data (default 0.05)
    :param isotope: string for filtering output dataframe such that columns from a particular isotope type are included
    (default 'light')
    :param logfile: path to file to log

    :return: a pandas dataframe with no index, columns of peptide_modified_sequence, rt_start (seconds),
    rt_end (seconds), charge, mz [note that this may only be approximate for modified peptides], protein_IDs
    """
    log('reading peptides from skyline dataframe ' + sample, logfile)

    # generate pandas df, filter, parse peptides, save encyc_rts values, save charge values, calculate mz, save proteins
    if q_value > 0.00:
        filt_entries = skyline_sub_df.loc[
            skyline_sub_df[isotope + ' ' + sample + ' ' + defs.DETECT_Q_VALUE_FIELD] < q_value]
    else:
        filt_entries = skyline_sub_df

    all_columns = filt_entries.columns
    key_columns = all_columns[(all_columns.str.contains(defs.RT_START_FIELD) |
                               all_columns.str.contains(defs.RT_END_FIELD))]
    pep_num_all = filt_entries.shape[0]
    log('found ' + str(pep_num_all) + ' peptides.', logfile)

    filt_entries.dropna(axis=0, how='any', subset=key_columns, inplace=True)
    peptides = [rename_peptide(i) for i in filt_entries[defs.PEPTIDE_MOD_SEQ_FIELD]]

    pep_num_clean = len(peptides)
    assert pep_num_clean > 0, "After filtering there were no peptides, aborting. Check your skyline output file."
    log('processing ' + str(pep_num_clean) + ' peptides.  ' + str(pep_num_all - pep_num_clean) +
        ' other peptides were dropped due to NaN RTs', logfile)

    rts_start = [round(t, 2) for t in filt_entries[isotope + ' ' + sample + ' ' + defs.RT_START_FIELD]]
    rts_end = [round(t, 2) for t in filt_entries[isotope + ' ' + sample + ' ' + defs.RT_END_FIELD]]
    zs = [i for i in filt_entries[defs.PEPTIDE_CHARGE_FIELD]]
    prot_names = [i for i in filt_entries[defs.PROTEIN_GENE_FIELD]]
    mzs = [mz for mz in filt_entries[isotope + ' ' + defs.ISOTOPE_FIND_FIELD]]
    pep_start = [map_int(start_res) for start_res in filt_entries[defs.START_POS_FIELD]]
    pep_end = [map_int(end_res) for end_res in filt_entries[defs.END_POS_FIELD]]

    # calc number of used peptides, and make sure all the arrays are the same size
    used_peptides = len(peptides)
    assert used_peptides == len(rts_start)
    assert used_peptides == len(rts_end)
    assert used_peptides == len(zs)
    assert used_peptides == len(prot_names)
    assert used_peptides == len(mzs)
    assert used_peptides == len(pep_start)
    assert used_peptides == len(pep_end)

    for rt_index in range(len(rts_start)):
        assert rts_start[rt_index] < rts_end[rt_index], log(
            'There is a problem with the retention times for peptide: ' + peptides[
                rt_index] + '. It has the following fields as starting and ending retention times: ' + str(
                rts_start[rt_index]) + ', ' + str(rts_end[rt_index]), logfile)

    # make a dataframe with your newly calculated entries
    pandas_dataframe = pd.DataFrame({peptides[i]: {'peptide_modified_sequence': peptides[i], 'rt_start': rts_start[i],
                                                   'rt_end': rts_end[i],
                                                   'charge': zs[i], 'mz': mzs[i], 'protein_IDs': prot_names[i],
                                                   'start_pos': pep_start[i], 'end_pos': pep_end[i]}
                                     for i in range(used_peptides)}).T
    pandas_dataframe = pandas_dataframe.sort_values(by=['rt_start'])
    pandas_dataframe = pandas_dataframe.dropna(axis=0, how='any', subset=['rt_start', 'rt_end', 'charge'])
    log('returning ' + str(pandas_dataframe.shape[0]) + ' peptides', logfile)
    return pandas_dataframe


def parse_skyline(path_to_skyline_csv, output_directory, sample_list=None, protein_list=None, isotope='light',
                  q_value=0.00, logfile=None):
    """
    Wrapper function to execute extract_skyline_sub and parse_sub_skyline on all proteins of interest. Useful if there
    are multiple samples in a given report file as it will then generate separate pysodist .tsv files for each of
    these samples.

    :param path_to_skyline_csv: string of directory path to location of skyline report
    :param output_directory: string with the directory where the output .csv files should be saved
    :param sample_list:  list of strings of individual sample identifiers in skyline report (default None results in
    analysis of all samples)
    :param protein_list: string of protein_gene_field identifier for particular proteins of interest (default None
    includes all proteins)
    :param isotope: string for filtering output dataframe such that columns from a particular isotope type are included
    (default 'light')
    :param q_value: float with a optional q_value cutoff to filter your data (default: no filtering)
    :param logfile: path to file to log

    :return: a list of pandas dataframes with no indices, columns of peptide_modified_sequence, rt_start (seconds),
    rt_end (seconds), charge,
    mz of the light species, protein_IDs, peptide_start_position, peptide_end_position
    """

    output_directory = output_directory.replace('\\', '/')
    if output_directory[-1] != '/':
        output_directory += '/'
    if not (os.path.exists(output_directory)):
        os.mkdir(output_directory)
    skyline_complete = pd.read_csv(path_to_skyline_csv, sep=',')
    if sample_list is None:
        sample_fields = [' '.join(i.split(defs.SAMPLE_FIND_FIELD)[0].split(' ')[1:]) for i in skyline_complete.columns
                         if defs.SAMPLE_FIND_FIELD in i]
        sample_list = [i.strip() for i in set(sample_fields)]
    output_list = []
    for sample in sample_list:
        log('working on sample: ' + sample, logfile)
        extracted_sub = extract_skyline_sub(skyline_complete, sample, protein_list=protein_list, isotope=isotope)
        output_list.append(parse_sub_skyline(extracted_sub, sample, q_value=q_value, isotope=isotope, logfile=logfile))
        write_directory = output_directory + sample
        try:
            os.mkdir(write_directory)
        except OSError:
            print('++++ The output directory: ' + write_directory +
                  ' already exists, and files within it may be overwritten. continue? [y/n]')
            choice = input().lower()
            if not choice == 'y':
                log('exiting parse_input as the desired output directory exists', logfile)
                raise
        log('Writing result to: ' + write_directory + '/pd_parsed_report.tsv')
        output_list[-1].to_csv(write_directory + '/pd_parsed_report.tsv', sep='\t',
                               columns=['rt_start', 'rt_end', 'peptide_modified_sequence', 'charge',
                                        'mz', 'protein_IDs', 'start_pos', 'end_pos'], index=False)
    return output_list


def add_args(parser):
    parser.add_argument('input', help='input file to parse.')
    parser.add_argument('--output_directory', default='./',
                        help='Output files will be saved in this folder: 1 directory '
                             'per sample in the skyline report. Default = ./')
    parser.add_argument('--sample_list', nargs='*', default=None,
                        help='An optional list of samples to parse. By default '
                             'all samples in the report are analyzed. Each sample separated by a space')
    parser.add_argument('--protein_list', nargs='*', default=None,
                        help='An optional list of the proteins to parse. By default, all proteins in the report are '
                             'analyzed. Each Protein Gene Name separated by a space.')
    parser.add_argument('--isotope', default='light', help='Be default, it is assumed that the report contains a light '
                                                           'isotope (no special labeling), if this field is not present'
                                                           'in the report, you can specify a different field here '
                                                           '(e.g. "heavy")')
    parser.add_argument('--q_value', default=0.00,
                        help='Used to optionally filter the report file based on the q_value. '
                             'By default, no q_value filtering is used.')
    parser.add_argument('--logfile', default=None, help='Optionally provide a path to a logfile to store outputs')
    return parser


def main(args):
    log('\n****INITIATING****', args.logfile)
    log('executed command: ' + " ".join(sys.argv), args.logfile)
    sample_list = args.sample_list
    parse_skyline(args.input, sample_list=sample_list, protein_list=args.protein_list,
                  isotope=args.isotope, q_value=args.q_value,
                  output_directory=args.output_directory, logfile=args.logfile)
    log('++++COMPLETED parse_input++++\n\n', args.logfile)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description='Pysodist input file parser. Used to parse report files from tools such'
                    ' as Skyline, EncyclopeDIA, or TPP, which provide a list of detected '
                    'peptides and retention times. Currently, only skyline report file '
                    'parsing is implemented.')
    add_args(argparser)
    main(argparser.parse_args())
