# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 15:29:11 2020

@author: joeyd
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

def write_scan(current_peptide, select_scan, file_name):
    '''
    Writes an individual scan to a file that can be read by default (fortran) isodist
    :param current_peptide: a pandas dataframe containg a single row of the parsed tsv file
    :param select_scan: a pandas dataframe with the selected scan (only including the desired m/z range). Can be interpolated, or summed, or raw.
    :param file_name: string pointing to the full path for the file to be written.

    returns: None
    '''
    select_scan.to_csv(file_name, sep=' ', header=False, index=False)
    return None

def extract_spectra(parsed_mzml, parsed_report, output_dir,
                    labeling='N15', save_interp_spectra=False, interp_res=0.001, sum_spectra_only=False, IO=True):
    '''
    Extracts spectra as specified in the parsed level 1 output. Writes individual spectra files as well as 
    a batch file for them. Optionally can interpolate the spectra, and sum them on a per-peptide basis.
    
    :param parsed_mzml: the parsed mzml file (should be output of parse_mzml)
    :param parsed_report: the full path to output file from the level1 parser. 
    :param output_dir: the full path to the 'base' directory. Within it, a spectra will be created.
    :param labeling: string defining the labeling type, which is used to determine the width of the m/z window. 
                     Currently implemented N15 and C13.
    :param save_interp_spectra: bool defining whether to save interpolated spectra instead of raw spectra (spectra are interpolated at interp_res)
    :param interp_res: the resolution to interpolate to when producing the summed spectra. Default 0.001.
    :param sum_spectra_only: boolean determining if only the summed spectra should be saved (single spectra per peptide). Uses the 
                        interp_res to do the summation. Default False.
    :param IO: boolean determining if IO should be printed

    returns: pandas dataframe ready to be written as a batch file(s).
    '''
    
    parsed_report_df = pd.read_csv(parsed_report,sep='\t')
    rt_array = parsed_mzml['retention_times']
    if IO:
        print('total spectra in mzml: '+ str(parsed_mzml['retention_times'].shape[0]))
        print('rt_range in mzml: ' + str(parsed_mzml['retention_times'][0]) + ' - ' + str(parsed_mzml['retention_times'][-1]) + ' mins')
        print('total peptides to extract: ' + str(parsed_report_df.shape[0]))
        print('peptide rt_range: ' + str(parsed_report_df['rt_start'].min()) + ' - ' + str(parsed_report_df['rt_end'].max()) + ' mins')
    spectra_dict = {}
    
    peaks_dir = output_dir+'spectra/'
    local_peaks_dir = './spectra/'
    try:
        os.mkdir(peaks_dir)
    except OSError:
        print('...the output spectra directory: ' + peaks_dir + ' already exists, and files within it may be overwritten. continue? [y/n]')
        choice = input().lower()
        if not choice=='y':
            raise
        
    for index, current_peptide in parsed_report_df.iterrows():
        
        first_scan = np.argmax(rt_array>=current_peptide['rt_start'])
        last_scan = np.argmax(rt_array>current_peptide['rt_end'])-1
        mz_range = mz_window(current_peptide['peptide_modified_sequence'], 
                             current_peptide['charge'], current_peptide['mz'], labeling=labeling)
        
        if IO:
            print('extracting ' + str(last_scan-first_scan) + ' spectra for peptide ' + str(index) + ' : ' + current_peptide['peptide_modified_sequence'])
        
        interp_mz_axis = np.arange(mz_range[0], mz_range[1], interp_res)
        interp_summed_intensity = np.zeros(interp_mz_axis.shape[0])
        
        peptide_base_name = "_".join([current_peptide['peptide_modified_sequence'],
                                      str(current_peptide['charge']),
                                      str(round(current_peptide['mz'],3)),
                                      str(current_peptide['start_pos'])+'-'+str(current_peptide['end_pos'])])
        
        for current_scan_num in range(first_scan, last_scan):
            scan_data = parsed_mzml['ms1_scans'][current_scan_num]
            select_scan_data = scan_data.loc[(scan_data['mz_data'] >= mz_range[0]) & (scan_data['mz_data'] <= mz_range[1])]
            
            interpolation = interp1d(scan_data['mz_data'], scan_data['intensity_data'])
            interp_intensity = interpolation(interp_mz_axis)
            interp_scan_data = pd.DataFrame({'mz_data':interp_mz_axis, 'intensity_data':interp_intensity})
            
            interp_summed_intensity += interp_intensity

            if sum_spectra_only is False:
                file_name = peptide_base_name+'_'+str(round(rt_array[current_scan_num],3))
                spectra_string = peaks_dir+file_name+'.tsv'
                local_spectra_string = local_peaks_dir+file_name+'.tsv'
                if save_interp_spectra:
                    write_scan(current_peptide, interp_scan_data, spectra_string)
                else:
                    write_scan(current_peptide, select_scan_data, spectra_string)
                    
                spectra_dict[file_name] = [current_peptide['peptide_modified_sequence'], 
                                         str(current_peptide['charge']),
                                         local_spectra_string]
                
        file_name =  peptide_base_name+'_SUM'
        spectra_string = peaks_dir+file_name+'.tsv'
        local_spectra_string = local_peaks_dir+file_name+'.tsv'
        interp_summed_spectra = pd.DataFrame({'mz_data':interp_mz_axis, 'intensity_data':interp_summed_intensity/((last_scan-first_scan)/2)})
        write_scan(current_peptide, interp_summed_spectra, spectra_string)
        spectra_dict[file_name]= [current_peptide['peptide_modified_sequence'], 
                                         str(current_peptide['charge']),
                                         local_spectra_string]
        
        spectra_df = pd.DataFrame.from_dict(spectra_dict, orient='index', columns=['peptide_modified_sequence', 'charge', 'spectra_file'])
        spectra_df.to_csv(output_dir+'pd_exported_peaks.tsv', sep='\t', index=False)
    return pd.DataFrame.from_dict(spectra_dict, orient='index', columns=['peptide_modified_sequence', 'charge', 'spectra_file'])

def parse_mzml(mzml_path, pickle_data=None):
    '''
    retrieves all scans from a portion of an mzml file and generates arrays of mz versus intensity at each RT

    :param mzml_path: string pointing to the .mzml file to extra scans from
    :param pickle_data: string pointing to a pickle file to save the extracted spectra into. Default is not saved

    :return: a dictionary with keys retention_times, and ms1_scans. retention_times points to a numpy array of each retention time.
    ms1_scans points to a list of pandas dataframes with columns mz_data and intensity_data. The index of this list corresponds to
    the position in the retention_times numpy array.
    '''
    with mzml.read(mzml_path) as mz_reader:
        scan_list=[]
        rt_list = []
        print('reading mzml file ' + mzml_path + '...')
        for scan in mz_reader:
            assert scan['ms level'] == 1
            mz_int_pd = pd.DataFrame({'mz_data':scan['m/z array'], 'intensity_data':scan['intensity array']})
            scan_list.append(mz_int_pd)
            new_rt = float(scan['scanList']['scan'][0]['scan start time'])
            rt_list.append(new_rt)   
    parsed_mz_file = {'retention_times': np.array(rt_list), 'ms1_scans': scan_list}
    if not(pickle_data is None):
        pickle.dump(parsed_mz_file, open(pickle_data, 'wb'))
    return parsed_mz_file

def mz_window(peptide_modified_sequence, peptide_charge, peptide_mz, labeling = 'N15', topomers_left=2, tomopers_right=5):
    '''
    calculates the mz window required to accomodate the range of mzs possible for all isotopes present in a given peptidde

    :param peptide_modified_sequence: string of peptide amino acid sequence including modifications
    :param peptide_charge: int of peptide charge state
    :param peptide_mz: float of mz for the unlabeled monoisotopic peptide
    :param labeling: isotope labeleing method ie. "K8R10" or "N15" (default None: unlabeled)

    :return: list of lower and upper range of possible mz values +/- 8 mass units for each peptide
    '''
    peptide_simple_seq = ''.join(re.split('\[\+\d+\.\d+\]', peptide_modified_sequence))
    if labeling == 'N15':
        n14_offset = defs.N15MASS-defs.N14MASS
        labeled_mz = ((sum([defs.AANITROGENS[i]*n14_offset for i in peptide_simple_seq]))/peptide_charge)+peptide_mz
    elif labeling == 'C13':
        c13_offset = defs.C13MASS-defs.C12MASS
        labeled_mz = ((sum([defs.AACARBONS[i]*c13_offset for i in peptide_simple_seq]))/peptide_charge)+peptide_mz    
    #Additional labeling methods may be entered here as needed

    else:
        labeled_mz = peptide_mz
    return [peptide_mz-(topomers_left/peptide_charge),labeled_mz+(tomopers_right/peptide_charge)]

def add_args(parser):
    parser.add_argument('mzml', help='the relative path to mzml file to be analyzed.')
    parser.add_argument('parsed_report', help='the parsed report file generated by pysodist parse_input.py')
    parser.add_argument('--labeling', default='N15', help='The labeling scheme used for the highest mass isotope envelope you expect to fit. E.g. N15 or C13')
    parser.add_argument('--interp_only', action='store_const', const=True, default=False, help='Only save the interpolated spectra instead of the raw spectra')
    parser.add_argument('--sum_only', action='store_const', const=True, default=False, help='Only save summed (and interpolated) spectra instead of all individual spectra. Results in 1 spectra per peptide')
    parser.add_argument('--interp_res', default=0.001, type=float, help='Set the interpolation delta m/z - typical values from 0.01 to 0.001')
    return parser

def main(args):
    parsed_report=args.parsed_report.replace('\\','/')
    output_dir = '/'.join(parsed_report.split('/')[:-1])+'/'
    assert(os.path.exists(parsed_report) is True)
    
    parsed_mzml = parse_mzml(args.mzml)
    extract_spectra(parsed_mzml, args.parsed_report, output_dir,
                               labeling=args.labeling, save_interp_spectra=args.interp_only,
                               interp_res = args.interp_res, sum_spectra_only=args.sum_only)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Pysodist .mzml parser, used to extract the relevant spectra for a single sample as defined \
                                     by the parsed report file (e.g. from Skyline or Encyclopedia) from the matching .mzml file. \
                                    For Thermo instruments, one should generate the .mzml file from the original .raw file using msconvert as follows: \
                                    .\msconvert.exe ".\[YOUR_RAW_FILE].raw" -o "./" --mzML --64 -v mz64 --inten32 --noindex --filter "msLevel 1" --zlib \
                                    This tool will extra the individual spectra as .tsv files (saved in a ./OUTPUT_peaks directory), \
                                        and product .batch and .in files to be used by Fortran isodist')
    add_args(parser)
    main(parser.parse_args())
