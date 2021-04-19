# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
"""
import warnings

warnings.filterwarnings("ignore")
import os
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from os import path
import argparse

warnings.filterwarnings("default")
import pysodist
import shutil
import pysodist.utils.skyline_report_defs as defs
from pysodist import utilities
log = utilities.log
vlog = utilities.vlog


def parse_isodist_csv(file_path):
    """
    Parse an isodist csv output file, and save it as a pandas dataframe

    :param file_path: the full path to the output csv file

    returns: a pandas dataframe with appropriate fields updated.
    """
    parsed_id_output = pd.read_csv(file_path, sep=',')
    parsed_name_fields = parsed_id_output['file'].str.split('_', expand=True)
    rt_column_index = parsed_name_fields.shape[1] - 1
    parsed_id_output['retention_time'] = parsed_name_fields[rt_column_index]
    parsed_id_output['CID'] = parsed_id_output['pep'] + '+' + parsed_id_output['z_charge'].astype(str)
    parsed_id_output['UID'] = parsed_id_output['CID'] + parsed_id_output['retention_time']
    all_CIDs = list(set(parsed_id_output['CID'].values))
    parsed_id_output['chisq'] = parsed_id_output['chisq'].astype(float)
    parsed_id_output['chisq'] = np.log10(parsed_id_output['chisq'])
    for CID in all_CIDs:
        all_RTs = np.array(
            [float(i) for i in parsed_id_output[parsed_id_output['CID'] == CID]['retention_time'].values if i != 'SUM'])
        min_RT = float(all_RTs.min())
        max_RT = float(all_RTs.max())
        range_RT = max_RT - min_RT
        parsed_id_output.loc[parsed_id_output['CID'] == CID, 'clean_RT'] = \
            parsed_id_output[parsed_id_output['CID'] == CID]['retention_time'].str.replace('SUM', str(max_RT + 0.01))
        parsed_id_output.loc[parsed_id_output['CID'] == CID, 'peak_position'] = \
            (parsed_id_output[parsed_id_output['CID'] == CID]['clean_RT'].astype(float) - min_RT) / range_RT

    return parsed_id_output


def read_spectra(isodist_result, isodist_fit_folder, working_path='./', extension='.tsv'):
    """
    Read a given fit spectra produced by isodist

    Parameters
    ----------
    isodist_result : single row of the isodist pandas dataframe corresponding to a single spectra
    isodist_fit_folder : string with the path to the fit folder
    working_path : string with the working path
    extension : string with the extension used in the original input spectra (result of extract spectra). Default=.tsv
    Optional, default of -1 does not use any interpolation

    Returns
    -------
    dict
        Returns a dictionary with the fields 'data':a pandas dictionary of the spectra;
        'name'; the name to put on the plots; 'isodist_output':the isodist dataframe row.

    """

    input_spectra = isodist_result['file']
    input_spectra = input_spectra.replace('./', working_path)
    fit_spectra = isodist_fit_folder + input_spectra.split('/')[-1] + '.fit'
    input_spectra = input_spectra + extension
    spectra_pd = pd.read_csv(input_spectra, sep=' ', header=None).rename(columns={0: 'm/z', 1: 'intensity'})
    fit_pd = pd.read_csv(fit_spectra, sep=',', header=None).rename(columns={0: 'mass', 1: 'intensity'})
    fit_pd['m/z'] = fit_pd['mass'] / isodist_result['z_charge']
    spectra_pd['fit_int'] = np.interp(spectra_pd['m/z'], fit_pd['m/z'], fit_pd['intensity'])
    spectra_pd['resid'] = spectra_pd['intensity'] - spectra_pd['fit_int']
    plottable_name = ' | '.join([isodist_result['protein'], isodist_result['pep'],
                                 '+' + str(isodist_result['z_charge']), isodist_result['retention_time'] + ' mins'])
    return {'data': spectra_pd, 'name': plottable_name, 'isodist_output': isodist_result}


def plot_fit(spectral_dict, main_axis, resid_axis, numerator=(defs.AMPU,), denominator=(defs.AMPU, defs.AMPF),
             ignore_fields=('file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID', 'clean_RT',
                            'peak_position', 'current_ratio'), fontsize=8, output=None):
    """
    Plot individual isodist fit on a dual axis, with the raw and fit data on the top, and the residual on the bottom.

    Parameters
    ----------
    spectral_dict : a dictionary with keys 'data','name', 'isodist_output'; typically returned from read_spectra.
    main_axis : a matplotlib axis to plot the primary data and the fit
    resid_axis : a matplotlib axis to plot the residuals
    numerator : a list of fields to put in the numerator for the calculated ratio (defaults to AMP_U), optional
    denominator : a list of fields to put in the denominator for the calculated ratio (defaults to [AMP_U, AMP_F),
    ignore_fields : a list of fields to exclude from the info box, optional.
    The default is ['file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID', 'clean_RT'].
    fontsize : int with the font size. The default is 8.
    output : a string describing the full path to the file to be saved, optional. The default is None.

    Returns
    -------
    A list of the two plotted axes (main_axis, and resid_axis)

    """
    sns.lineplot(data=spectral_dict['data'], x='m/z', y='fit_int', ax=main_axis, color='orange').set_title(
        spectral_dict['name'])
    sns.scatterplot(data=spectral_dict['data'], x='m/z', y='intensity', marker='x', ax=main_axis, color='black')
    sns.lineplot(data=spectral_dict['data'], x='m/z', y='resid', ax=resid_axis, color='red')
    result_string = _parse_result_string(spectral_dict['isodist_output'], numerator=numerator, denominator=denominator,
                                         ignore_fields=ignore_fields)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    main_axis.text(0.2, 0.95, result_string, transform=main_axis.transAxes, fontsize=fontsize, verticalalignment='top',
                   bbox=props)

    if output is not None:
        plt.savefig(output)
        plt.close()

    return [main_axis, resid_axis]


def _parse_result_string(isodist_row, numerator=(defs.AMPU,), denominator=(defs.AMPU, defs.AMPF),
                         ignore_fields=('file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID',
                                        'clean_RT', 'peak_position', 'current_ratio')):
    """
    Helper function to generate a result string for plotted spectra/fits


    Parameters
    ----------
    isodist_row : the row from an isodist pandas dataframe (result of parse_isodist_csv)
    numerator : a list of fields to put in the numerator for the calculated ratio (defaults to AMP_U), optional
    denominator : a list of fields to put in the denominator for the calculated ratio (defaults to [AMP_U, AMP_F)
    ignore_fields : a list of fields to exclude from the info box, optional.
    The default is ['file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'CID', 'clean_RT', 'UID'].

    Returns
    -------
    A string with the summarized isodist results for this fit to be placed on the plot.

    """

    ratio = str(round(np.array([isodist_row[num] for num in numerator]).sum() / np.array(
        [isodist_row[den] for den in denominator]).sum(), 3))
    all_fields = list(isodist_row.index.values)
    fields = [field for field in all_fields if field not in ignore_fields]
    output_string = '\n'.join([field + ": " + _clean_round(isodist_row[field], 3) for field in fields])
    output_string = output_string + '\n[' + '+'.join([num for num in numerator]) + ']/[' + '+'.join(
        [den for den in denominator]) + '] = ' + ratio

    return output_string


def _clean_round(input_value, decimal_places):
    """
    Rounds an input value to the desired number of decimal places; handles strings and floats

    Parameters
    ----------
    input_value : string or float of the number to be rounded.
    decimal_places : int of the number of decimal places

    Returns
    -------
    a string with the rounded value.

    """
    try:
        cleaned_input = str(round(float(input_value), decimal_places))
    except ValueError:
        cleaned_input = str(input_value)
    return cleaned_input


def parse_raw(batch_row):
    protein, _, mz, rt = batch_row['spectra'].split('/')[-1].split('_')
    rt = rt.split('.txt')[0]
    fit_file = batch_row['spectra'].split('.txt')[0] + '.fit'
    dat_file = batch_row['spectra'].split('.txt')[0] + '.dat'
    assert path.exists(fit_file)
    assert path.exists(dat_file)
    return {'protein_name': protein, 'mz': float(mz), 'retention_time': float(rt),
            'raw_data': batch_row['spectra'], 'fit_data': fit_file, 'dat_file': dat_file,
            'charge': batch_row['charge'], 'pep_seq': batch_row['pep_seq']}


def plot_single_spectra(id_row, fit_folder, working_path='./', extension='.tsv', fig_size=(4, 6),
                        fontsize=6,
                        numerator=(defs.AMPU,), denominator=(defs.AMPU, defs.AMPF),
                        ignore_fields=('file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID',
                                       'clean_RT', 'peak_position', 'current_ratio')):
    """
    Plots a single spectra from a given isodist dataframe row (produces the full figure and axes)

    Parameters
    ----------
    id_row : the row from an isodist pandas dataframe (result of parse_isodist_csv)
    fit_folder : string of the full path to the folder containing the fit files.
    working_path : string with the full path to the working directory. Optional, default = './'
    extension : string with the spectra file type extension (default is .tsv).
    fig_size : tuple of floats with the zsize of the resulting figures. The default is (4,6).
    fontsize : float with the font size. The default is 6.
    numerator : a list of fields to put in the numerator for the calculated ratio (defaults to AMP_U), optional
    denominator : a list of fields to put in the denominator for the calculated ratio (defaults to [AMP_U, AMP_F)
    ignore_fields : a list of fields to exclude from the info box, optional.
    The default is ['file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID', 'clean_RT'].

    Returns
    -------
    A list with the main figure axis and the residual figure axis.

    """
    spectral_dict = read_spectra(id_row, fit_folder, working_path=working_path, extension=extension)
    fig, axes = plt.subplots(ncols=1, nrows=2, gridspec_kw={'height_ratios': [5, 1]}, figsize=fig_size, sharex=True)
    main_plot, resid_plot = plot_fit(spectral_dict, axes[0], axes[1], numerator=numerator, denominator=denominator,
                                     ignore_fields=ignore_fields, fontsize=fontsize)

    return [main_plot, resid_plot]


def plot_spectra_group(id_output, fit_folder, working_path='./', page_size=(22, 14), fig_layout=(3, 3),
                       numerator=(defs.AMPU,), denominator=(defs.AMPU, defs.AMPF),
                       ignore_fields=('file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID',
                                      'clean_RT', 'peak_position', 'current_ratio'),
                       saved_output=None, png=True, pdf=False, logfile=None):
    """
    Plots a series of individual spectra and fits. Typically these would all be from a related peptide or proteins,
    but in theory they can be any set of spectra as defined by rows of the isodist pandas dataframe.

    Parameters
    ----------
    id_output : selected rows from the parsed isodist dataframe corresponding to the spectra to be plotted.
    fit_folder : string of the full path to the folder containing the fit files.
    working_path : string with the working path.
    page_size : a tuple with the size of each page to be generated. Optional with a default of (22,14)
    numerator : a list of fields to put in the numerator for the calculated ratio (defaults to AMP_U), optional
    denominator : a list of fields to put in the denominator for the calculated ratio (defaults to [AMP_U, AMP_F)
    ignore_fields : a list of fields to exclude from the info box, optional.
    The default is ['file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID', 'clean_RT'].
    fig_layout : tuple with the layout for the panels per page. Optional with a default of (3,3)
    saved_output : string with the full path to an example  file to be generated (without the desired extension).
    Each will be created as this name_page#_. Optional, default is None, and no files will be saved.
    png : bool defining if .png files are saved. Optional (default is True)
    pdf : bool defining if .pdf files are saved. Optional (default is False)
    logfile : path to file to log

    Returns
    -------
    None.

    """

    indices = id_output.index
    num_spectra = indices.shape[0]
    fig_per_page = fig_layout[0] * fig_layout[1]
    num_pages = num_spectra // fig_per_page + 1

    fig_cols = fig_layout[0]
    fig_rows = fig_layout[1]

    for page in range(num_pages):
        log('plotting page ' + str(page) + ' of ' + str(num_pages - 1) + '...', logfile)
        fig, axes = plt.subplots(ncols=fig_cols, nrows=fig_rows * 2,
                                 gridspec_kw={'height_ratios': [5, 1] * fig_layout[1]}, figsize=page_size)
        for counter, row in enumerate(indices[page * fig_per_page:(page + 1) * fig_per_page]):
            plottable_dict = read_spectra(id_output.loc[row], fit_folder, working_path=working_path)

            main_col = (counter % fig_cols)
            resid_col = main_col

            main_row = (counter // fig_cols) * 2
            resid_row = main_row + 1

            plot_fit(plottable_dict, axes[main_row, main_col], axes[resid_row, resid_col], numerator=numerator,
                     denominator=denominator, ignore_fields=ignore_fields)
        plt.tight_layout()
        if saved_output is not None:
            current_page_output = saved_output + '_' + str(page)
            if png:
                plt.savefig(current_page_output + '.png')
            if pdf:
                plt.savefig(current_page_output + '.pdf')
            plt.close()


def get_by_peptide(id_output, peptide_sequence, charge=None):
    """
    Get all associated spectra given a peptide sequence.

    Parameters
    ----------
    id_output : full isodist pandas dataframe pandas dataframe (result of parse_isodist_csv)
    peptide_sequence : String of the peptide to use to find all related peptides.
    charge : int with the optional charge state to search for.
    Optional. If None provided, all spectra of all charges states returned.

    Returns
    -------
    peptide_subset : a pandas dataframe with the subset of related spectral rows.

    """
    peptide_subset = id_output[id_output['UID'].str.contains(peptide_sequence)]
    if not (charge is None):
        peptide_subset = peptide_subset[peptide_subset['UID'].str.contains(str(charge) + '_')]
    return peptide_subset


def get_by_protein(id_output, protein_name):
    """
    Helper function to get the rows containing the protein_name.

    Parameters
    ----------
    id_output : full isodist pandas dataframe pandas dataframe (result of parse_isodist_csv)
    protein_name : string with the name of the protein to get by

    Returns
    -------
    the subset of the input dataframe that contains the protein_name in the 'protein field'

    """
    # to get just the list of peptides: list(set(get_by_protein(id_output, 'ENO2_YEAST')['pep'].values))
    return id_output[id_output['protein'].str.contains(protein_name)]


def plot_ratios(related_spectra, y_label, marker_size=6, palette='winter', fig_height=4):
    """
    Plots a single axis with a series of related spectra - typically all of the spectra for a given protein.
    They are broken into groups based on the peptide field, and the y-value is whatever is in "current_ratio";
    this should be calculated BEFORE using this function. Dots are separated along the 'peak_position' field.

    Parameters
    ----------
    related_spectra :  a pandas dataframe with the subset of related spectral rows.
    y_label : a string describing what the current ratio field is reporting.
    marker_size : float, optional. The default is 6.
    palette : string, optional string for the color pallet to use along the peak position. The default is 'winter'.
    fig_height : float wit with the figure height to plot. Default=4.

    Returns
    -------
    figure : the resulting matplotlib figure

    """

    all_peptides = list(set(related_spectra['CID'].values))
    all_peptides.sort()
    num_peps = len(all_peptides)
    fig, axes = plt.subplots(nrows=1, ncols=num_peps, figsize=[fig_height * (2 * num_peps / 3), fig_height],
                             sharey=True, sharex=True)

    for count, peptide in enumerate(all_peptides):
        if num_peps > 1:
            current_axis = axes[count]
        else:
            current_axis = axes
        peptides_df = related_spectra[related_spectra['CID'] == peptide]
        current_axis = sns.scatterplot(data=peptides_df, x='peak_position', y='current_ratio', hue='peak_position',
                                       ax=current_axis, palette=palette, markersize=marker_size)
        current_axis.set_xlabel(peptide)
        current_axis.get_legend().remove()
        current_axis.set_ylabel(y_label)
    # sns.swarmplot(x='pep', y='current_ratio', hue='peak_position', data=related_spectra,
    #              dodge=True, size=marker_size, ax=axis, palette=palette) #TOO SLOW
    # sns.stripplot(x='pep', y='current_ratio', hue='peak_position', data=related_spectra,
    #              dodge=True, size=marker_size, ax=axis, palette=palette) #TOO SLOW
    # sns.boxplot(x='pep', y='current_ratio', data=related_spectra, ax=axis, color='w')
    # sns.violinplot(x='pep', y='current_ratio', data=related_spectra, innter=None, ax=axis)
    return fig


def plot_all_ratios(id_output, numerator=(defs.AMPU,), denominator=(defs.AMPU, defs.AMPF), fig_height=4,
                    saved_output_path=None, png=True, pdf=False, protein_list=None, logfile=None):
    """
    Plots the observed ratios for each spectra for each protein in the input dataframe.

    Parameters
    ----------
    id_output : full isodist pandas dataframe pandas dataframe (result of parse_isodist_csv)
    numerator : a list of fields to put in the numerator for the calculated ratio (defaults to AMP_U), optional
    denominator : a list of fields to put in the denominator for the calculated ratio (defaults to [AMP_U, AMP_F)
    fig_height : float with the height of each field. Optional, default is 4.
    saved_output_path : string with path to save figures (will be path+protein_name+extension).
    Optional, default None is to not save the figures
    png : bool defining if .png files are saved. Optional (default is True)
    pdf : bool defining if .pd
    logfile : path to file where log files are saved. Optional (default is False)

    protein_list : list of strings of the proteins to plot.
    Optional, default is None, which results in plotting for all proteins in the id_output dataframe.

    Returns
    -------
    None.

    """

    y_label = '[' + '+'.join([num for num in numerator]) + ']/[' + '+'.join([den for den in denominator]) + ']'
    id_output = set_current_ratio(id_output, numerator=numerator, denominator=denominator)
    if protein_list is None:
        protein_list = list(set(id_output['protein'].values))
    protein_list.sort()
    for protein in protein_list:
        log('plotting abundance for protein: ' + protein, logfile)
        related_spectra = get_by_protein(id_output, protein)
        fig = plot_ratios(related_spectra, y_label, marker_size=fig_height + 1, fig_height=fig_height)
        fig.suptitle(protein)
        plt.tight_layout()
        if saved_output_path is not None:
            save_name = saved_output_path + protein + '_abundance'
            if png:
                plt.savefig(save_name + '.png')
            if pdf:
                plt.savefig(save_name + '.pdf')
            plt.close()


def set_current_ratio(id_output, numerator=(defs.AMPU,), denominator=(defs.AMPU, defs.AMPF)):
    """
    Calculates an abundance ratio, returns the updated pandas dataframe

    Parameters
    ----------
    id_output : full isodist pandas dataframe pandas dataframe (result of parse_isodist_csv)
    numerator : a list of fields to put in the numerator for the calculated ratio (defaults to [AMP_U]), optional
    denominator : a list of fields to put in the denominator for the calculated ratio (defaults to [AMP_F), optional

    Returns
    -------
    the resulting pandas dataframe.

    """
    id_output.loc[:, 'current_ratio'] = (id_output[numerator].sum(axis=1) / id_output[denominator].sum(axis=1)).round(3)
    return id_output


def plot_csv_stats(id_output, output_path=None, png=True, pdf=False,
                   histo_fields=('chisq', 'B', 'OFF', 'GW', 'current_ratio', 'FRC_NX', 'FRC_NY'), figsize=(10, 12),
                   bins=100):
    """
    Plots histograms of fields from the parsed pandas dataframe.
    Also plots histograms for numerator/denominator, and numerator/(numerator+denominator) for each spectra.

    Parameters
    ----------
    id_output : full isodist pandas dataframe pandas dataframe (result of parse_isodist_csv)
    output_path : string of the full path to the file to save (should include extension). Optional, default is None
    png : bool to output a .png file (default = False)
    pdf : bool to output a .pdf file (default = False)
    histo_fields : list of strings with the isodist fields to plot histograms for.
    Optional, default is ['chisq', 'B', 'OFF', 'GW', 'FRC_NX', 'FRC_NY'].
    figsize : tuple with the desired figure size. Optional, default is [10,12].
    bins : int for the number of bins in the histograms. Optional, default is 100.

    Returns
    -------
    the resulting figure.

    """
    histo_fields = [h for h in histo_fields if h in id_output.columns]
    num_rows = len(histo_fields)
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=figsize)
    for row, field in enumerate(histo_fields):
        sns.histplot(id_output, x=field, ax=axes[row, 0], bins=bins)
        sns.histplot(id_output, x=field, ax=axes[row, 1], bins=bins)
        median = id_output[field].median()
        std = id_output[field].std()
        axes[row, 1].set_xlim([median - std, median + std])

    plt.tight_layout()
    if not (output_path is None):
        if png:
            plt.savefig(output_path + 'summary.png')
        if pdf:
            plt.savefig(output_path + 'summary.pdf')
        plt.close()
    return fig


def add_args(parser):
    parser.add_argument('input_file', help='path to the compiled isodist .csv file with all results (1 row/fit).')
    parser.add_argument('fit_folder', help='path to the folder containing all of the isodist .fit files.')
    parser.add_argument('output_folder',
                        help='path to a folder to save all of the fits into. Will be created if it does not exist.')
    parser.add_argument('--numerator', nargs='+', default=['AMP_U'],
                        help='list of the fields to use in the numerator of the abundance ratio calculation '
                             '(typically AMP_U, AMP_L, AMP_F, or some combination. Default is AMP_U')
    parser.add_argument('--denominator', nargs='+', default=['AMP_U', 'AMP_F'],
                        help='list of the fields to use in the denominator of the abundance ratio calculation '
                             '(typically AMP_U, AMP_L, AMP_F, or some combination. Default is AMP_U, AMP_F')
    parser.add_argument('--no_png', action='store_const', const=True, default=False,
                        help='By default .png files for the plots will be saved. '
                             'This option forces these to not be saved.')
    parser.add_argument('--no_pdf', action='store_const', const=True, default=False,
                        help='By default .pdf files for the plots will be saved. '
                             'This option forces these to not be saved.')
    parser.add_argument('--logfile', default=None, help='Optionally provide a path to a logfile to store outputs')

    return parser


def main(args):
    input_file = args.input_file.replace('\\', '/')
    working_path = '/'.join(input_file.split('/')[:-2]) + '/'
    fit_folder = args.fit_folder.replace('\\', '/')
    logfile = args.logfile.replace('\\', '/')
    if fit_folder[-1] != '/':
        fit_folder += '/'
    output_folder = args.output_folder.replace('\\', '/')
    if output_folder[-1] != '/':
        output_folder += '/'

    log('parsing isodist csv file: ' + input_file, logfile)
    isodist_output = parse_isodist_csv(input_file)

    all_isodist_columns = isodist_output.columns
    for num in args.numerator:
        assert num in all_isodist_columns, 'provided numerator ' + num + ' not present in report_file'
    for den in args.denominator:
        assert den in all_isodist_columns, 'provided denominator ' + den + ' not present in report_file'

    current_ratio_string = '[' + '+'.join([num for num in args.numerator]) + ']/[' + '+'.join(
        [den for den in args.denominator]) + ']'
    log('all of the following plots will use current ratio as: ' + current_ratio_string, logfile)
    isodist_output = set_current_ratio(isodist_output, numerator=args.numerator, denominator=args.denominator)

    try:
        os.mkdir(output_folder)
        # copy over template if file doesn't exist
    except OSError:
        print(
            '...the output directory: ' + output_folder +
            ' already exists, and files within it may be overwritten. continue? [y/n]')
        choice = input().lower()
        if not choice == 'y':
            raise

    log('creating jupyter notebook for interactive analysis...', logfile)
    source = f'{pysodist._ROOT}/utils/analysis_template.ipynb'
    shutil.copyfile(source, output_folder + '/analysis_template.ipynb')

    isodist_output.to_csv(output_folder + 'isodist_result.csv')

    assert not (args.no_png and args.no_pdf)

    log('using fit spectra from directory: ' + fit_folder, logfile)
    all_proteins = list(set(isodist_output['protein'].values))

    log('saving plots for ' + str(len(all_proteins)) + ' proteins to: ' + output_folder, logfile)
    all_proteins.sort()
    for protein in all_proteins:
        log('plotting spectra for protein: ' + protein, logfile)
        related_spectra = get_by_protein(isodist_output, protein)
        plot_spectra_group(related_spectra, fit_folder, working_path=working_path, numerator=args.numerator,
                           denominator=args.denominator, png=not args.no_png,
                           pdf=not args.no_pdf, saved_output=output_folder + protein, logfile=logfile)

    log('calculating and plotting abundances...', logfile)
    plot_all_ratios(isodist_output, numerator=args.numerator, denominator=args.denominator,
                    saved_output_path=output_folder, png=not args.no_png, pdf=not args.no_pdf, logfile=logfile)
    log('plotting the csv stat histograms...', logfile)
    plot_csv_stats(isodist_output, output_path=output_folder, png=not args.no_png,
                   pdf=not args.no_pdf)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description='Pysodist plotter to generate and save plots of the isodist fits. '
                                                    'Also saves datastructures to be analyzed using jupyterlab')
    add_args(argparser)
    main(argparser.parse_args())
