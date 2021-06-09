import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pysodist.commands.plot_spectra as ps
from scipy import optimize
from multiprocessing import Pool
from itertools import chain


def add_args(parser):
    parser.add_argument('--spectrassr', type=int, required=True, help='Spectra SSR to filter scans by')
    parser.add_argument('--ssr', type=int, required=True, help='Gaussian fitting SSR to filter peptides by')
    parser.add_argument('--threads', type=int, required=True, help='Number of threads to use in parallel')
    parser.add_argument('--workdir', type=str, default='./', help='Working directory')
    parser.add_argument('--fitdir', type=str, help='Directory in which isodist fits are stored')
    parser.add_argument('--outdir', type=str, help='Directory for storing output plots and files')
    parser.add_argument('--isocsv', type=str, help='Output .csv file from isodist')
    parser.add_argument('--expect', nargs='+', default=None, type=float, help='Expectred ["U/F", "L/F", "U/L"] ratios')
    parser.add_argument('--filt', type=str, default='AMP_F',
                        help='Value to fit and filter by (AMP_U, AMP_L, AMP_F, or sum). Default is AMP_F')
    parser.add_argument('--plotfilt', type=bool, default=True, help='Whether to plot retained and filtered-out spectra')
    return parser


def parse_isodist_csv(file_path):
    parsed_id_output = pd.read_csv(file_path, sep=',')
    drop_inds = [i for i in range(0, len(parsed_id_output)) if
                 (parsed_id_output['file'][i].count('.') < 3) & ('SUM' not in parsed_id_output['file'][i])]
    parsed_id_output = parsed_id_output.drop(drop_inds)
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
        if len(all_RTs) > 0:
            min_RT = float(all_RTs.min())
            max_RT = float(all_RTs.max())
            range_RT = max_RT - min_RT
            parsed_id_output.loc[parsed_id_output['CID'] == CID, 'clean_RT'] = \
            parsed_id_output[parsed_id_output['CID'] == CID]['retention_time'].str.replace('SUM', str(max_RT + 0.01))
            parsed_id_output.loc[parsed_id_output['CID'] == CID, 'peak_position'] = (parsed_id_output[parsed_id_output[
                                                                                                          'CID'] == CID][
                                                                                         'clean_RT'].astype(
                float) - min_RT) / range_RT

    return parsed_id_output


def piecewise_triangle(x, a, c, s, w):
    return np.piecewise(x, [abs(x + s) < w, abs(x + s) >= w], [lambda x: a - a * abs((x + s) / w) + c, lambda x: c])


def gaussian(x, a, c, s, w):
    return a * np.exp(-(x + s) ** 2 / (2 * w ** 2)) + c


def fit_amps(peptide_df, amp_col, func):
    x = pd.to_numeric(peptide_df['retention_time']).values

    if type(amp_col) == str:
        y = peptide_df[amp_col].values
    if type(amp_col) == list:
        y = peptide_df[amp_col[0]]
        for i in range(1, len(amp_col)):
            y = y + peptide_df[amp_col[i]]
        y = y.values

    if func == piecewise_triangle or func == gaussian:
        c_guess = y.min()
        if len(y) > 1:
            if max(y) <= 2 * np.sort(y)[-2]:
                s_guess = -1 * x[np.argmax(y)]
                a_guess = y.max() - y.min()
            elif max(y) > 2 * np.sort(y)[-2]:
                s_guess = -1 * x[np.where(y == np.sort(y)[-2])[0][0]]
                a_guess = np.sort(y)[-2] - y.min()
        else:
            s_guess = -1 * (x.min() + (x.max() - x.min()) / 2)
            a_guess = y.max() - y.min()
        w_guess = (x.max() - x.min()) / 3
        guesses = [a_guess, c_guess, s_guess, w_guess]
        bounds = ([0, -np.inf, -1 * (x.max()), 0], [np.inf, np.inf, -1 * (x.min()), x.max() - x.min()])

    popt, pcov = optimize.curve_fit(func, x, y, p0=guesses, bounds=bounds)
    return x, y, popt, pcov


def calc_ssr(x, y, func, popt):
    resid_y = func(x, *popt)
    y_mean = np.mean(y)
    ssr_adj = sum(((resid_y - y) / y_mean) ** 2) / len(resid_y)
    return ssr_adj


def plot_fit_mod(spectral_dict, main_axis, resid_axis=None, numerator=['AMP_U'], denominator=['AMP_U', 'AMP_F'],
                 info_box=True,
                 ignore_fields=['file', 'protein', 'pep', 'mw', 'z_charge', 'retention_time', 'UID', 'CID', 'clean_RT',
                                'peak_position', 'current_ratio'], fontsize=8, output=None):
    sns.lineplot(data=spectral_dict['data'], x='m/z', y='fit_int', ax=main_axis, color='orange').set_title(
        spectral_dict['name'])
    sns.scatterplot(data=spectral_dict['data'], x='m/z', y='intensity', marker='x', ax=main_axis, color='black')
    if resid_axis:
        sns.lineplot(data=spectral_dict['data'], x='m/z', y='resid', ax=resid_axis, color='red')
    if info_box:
        result_string = ps._parse_result_string(spectral_dict['isodist_output'], numerator=numerator,
                                                denominator=denominator, ignore_fields=ignore_fields)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        main_axis.text(0.2, 0.95, result_string, transform=main_axis.transAxes, fontsize=fontsize,
                       verticalalignment='top', bbox=props)

    if not output is None:
        plt.savefig(output)
        plt.close()

    return [main_axis, resid_axis]


def run_fit_amps(inputs):
    fails = []
    pep_df, fit_col, fit_func, fig_output, ft, wd = inputs
    pep = pep_df['pep'].unique()[0]
    popt = None
    adj_ssr = None

    try:
        nrows = int(np.ceil((len(pep_df) + 1) / 2))

        page_height = 4 * nrows
        fig, ax = plt.subplots(nrows, 2, figsize=(12, page_height), tight_layout=True)
        ax = ax.flatten()

        x_vals, y_vals, popt, pcov = fit_amps(pep_df, fit_col, fit_func)
        adj_ssr = calc_ssr(x_vals, y_vals, fit_func, popt)
        fitx = np.linspace(x_vals.min(), x_vals.max(), 100)
        fity = fit_func(fitx, *popt)

        ax[0].plot(fitx, fity, '--r')
        ax[0].scatter(x_vals, y_vals, s=10)

        ax[0].set_title(pep)

        for j in range(0, len(pep_df)):
            id_row = pep_df.iloc[j]
            spectral_dict = ps.read_spectra(id_row, ft, working_path=wd)
            main, resid = plot_fit_mod(spectral_dict, ax[j + 1], info_box=False)

        fig.savefig(fig_output + pep + '.png')
        plt.close()
    except TypeError:
        fails.append(pep)
    except RuntimeError:
        fails.append(pep)
    except ValueError:
        fails.append(pep)

    return pep, popt, adj_ssr, fails


def plot_fits_spectra(df, fit_col, fit_func, fit_output, fit_output_cols, fig_output, ft, num_threads, wd):
    all_peps = df['pep'].unique()
    pep_dfs_list = []

    for pep in all_peps:
        sub = df[df['pep'] == pep]
        pep_dfs_list.append([sub, fit_col, fit_func, fig_output, ft, wd])

    pool = Pool(num_threads)
    fits = pool.map(run_fit_amps, pep_dfs_list)

    for i in fits:
        pep, popt, adj_ssr = i[0:3]
        inds = df[df['pep'] == pep].index
        fit_output.loc[inds, fit_output_cols[1:]] = popt
        fit_output.loc[inds, fit_output_cols[0]] = adj_ssr

    fail_list = [i[3] for i in fits if len(i[3]) > 0]
    fail_list = list(chain(*fail_list))

    return fit_output, fail_list


def calc_spectra_ssr(inputs):
    df, row, ft, wd = inputs
    spectra = ps.read_spectra(df.loc[row], ft, working_path=wd)
    spectra_ssr = np.sum((spectra['data']['resid'] / np.mean(spectra['data']['intensity'])) ** 2) / len(
        spectra['data']['resid'])
    return row, spectra_ssr


def IQR_boxplot(df, expecteds=None, scatter_x=None, color='b', ax=None, label=None):
    if not scatter_x:
        scatter_x = np.linspace(0, len(df.columns) - 1, len(df.columns))
    scatter_y = []
    scatter_error = np.empty([2, len(df.columns)])

    i = 0
    for col in df.columns:
        q75, q50, q25 = np.percentile(df[col].dropna(), [75, 50, 25])
        scatter_y.append(q50)
        scatter_error[0, i] = q25
        scatter_error[1, i] = q75
        i = i + 1

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    fmt = color + 'o'
    ax.errorbar(scatter_x, scatter_y, yerr=scatter_error, fmt=fmt, capsize=4, label=label)

    if expecteds is not None:
        for i in range(0, len(expecteds)):
            ax.plot([scatter_x[i] - 0.5, scatter_x[i] + 0.5], [expecteds[i], expecteds[i]], '--r')

    return


def plot_filt_spectra(inputs):
    pep_df, output, ft, wd = inputs
    pep = pep_df['pep'].unique()[0]
    nrows = int(np.ceil(len(pep_df) / 2))
    fig, ax = plt.subplots(nrows, 2, figsize=(12, 4 * nrows))
    ax = ax.flatten()

    i = 0
    for spectra in pep_df.index:
        spectral_dict = ps.read_spectra(pep_df.loc[spectra], ft, working_path=wd)
        plot_fit_mod(spectral_dict, ax[i], info_box=False)
        i = i + 1

    fig.savefig(output + pep + '.png')
    plt.close()
    return


def main(args):
    print('loading data...')
    wd = args.workdir

    if args.fitdir is not None:
        ft = args.fitdir
    else:
        ft = wd + [i for i in os.listdir(wd) if '_isodist_fits' in i][0] + '/'

    if args.outdir is not None:
        out = args.outdir
    else:
        out = wd + 'gaussian_' + str(args.spectrassr) + '_' + str(args.ssr) + '/'
    if not os.path.exists(out):
        os.mkdir(out)

    if args.isocsv is not None:
        isodist_result = args.isocsv
    else:
        fitname = [i for i in os.listdir(wd) if '_isodist_outputs' in i][0]
        isodist_result = wd + fitname + '/' + fitname.split('_isodist')[0] + '_output.csv'

    id_result = parse_isodist_csv(isodist_result)
    id_result.loc[:, 'U/F'] = (id_result['AMP_U'] / id_result['AMP_F']).round(3)
    id_result.loc[:, 'L/F'] = (id_result['AMP_L'] / id_result['AMP_F']).round(3)
    id_result.loc[:, 'U/L'] = (id_result['AMP_U'] / id_result['AMP_L']).round(3)

    ratio_cols = ['U/F', 'L/F', 'U/L']

    nonsum_id = id_result[id_result['retention_time'] != 'SUM']
    print('Calculating spectra SSR...')

    spectra_ssr_list = []
    for i in nonsum_id.index:
        spectra_ssr_list.append([nonsum_id, i, ft, wd])

    pool = Pool(args.threads)
    spectra_ssr_results = pool.map(calc_spectra_ssr, spectra_ssr_list)

    for i in spectra_ssr_results:
        row, spectra_ssr = i
        nonsum_id.loc[row, ['Spectra SSR']] = spectra_ssr
    nonsum_id.to_csv(out + 'nonsum_id.csv')

    print('Filtering by GW and amplitude...')
    gw_filt = nonsum_id[(nonsum_id['GW'] >= 0.01) & (nonsum_id['GW'] <= 1)]
    gw_filt = gw_filt[gw_filt['Spectra SSR'] <= np.percentile(gw_filt['Spectra SSR'], args.spectrassr)]
    spectra = []
    all_peps = gw_filt['pep'].unique()
    for pep in all_peps:
        pep_df = gw_filt[gw_filt['pep'] == pep]
        pep_meds = pep_df.median()
        sub_pep = pep_df[(pep_df['AMP_U'] >= 0.1 * pep_meds['AMP_U']) & (pep_df['AMP_L'] >= 0.1 * pep_meds['AMP_L']) & (
                    pep_df['AMP_F'] >= 0.1 * pep_meds['AMP_F']) & (pep_df['AMP_U'] <= 20 * pep_meds['AMP_U']) & (
                                     pep_df['AMP_L'] <= 20 * pep_meds['AMP_L']) & (
                                     pep_df['AMP_F'] <= 20 * pep_meds['AMP_F'])]
        if len(sub_pep) > 3:
            for i in range(0, len(sub_pep.index)):
                spectra.append(sub_pep.index.to_list()[i])

    gw_amp_filt = gw_filt.loc[spectra]
    gw_amp_filt.to_csv(out + 'gw_amp_filt.csv')

    fit_results = pd.DataFrame(index=gw_amp_filt.index)

    abbrev = str(args.filt).split('_')[1].lower()
    result_cols = ['SSR_gaussian_' + abbrev, 'Gaussian_' + abbrev + '_a', 'Gaussian' + abbrev + '_c',
                   'Gaussian_' + abbrev + '_s', 'Gaussian_' + abbrev + '_w']
    storefits = out + 'Gaussian_' + abbrev + '/'
    os.mkdir(storefits)
    print('fitting and plotting ' + str(args.filt) + '...')
    fit_results, failed_fits = plot_fits_spectra(gw_amp_filt, args.filt, gaussian, fit_results, result_cols, storefits,
                                                 ft, args.threads, wd)

    print('writing fit_results file...')
    fit_results.to_csv(out + 'fit_results.csv')

    print('writing failed fits...')
    fail_file = out + 'failed_fits.txt'
    with open(fail_file, 'w') as f:
        f.write('Failed fits\n')
        for j in failed_fits:
            f.write('\t' + str(j) + '\n')

    fit_results = fit_results.fillna(np.inf)
    f_filt = fit_results[
        fit_results['SSR_gaussian_' + abbrev] <= np.percentile(fit_results['SSR_gaussian_' + abbrev], args.ssr)].index

    expected_ratios = args.expect
    subf = gw_amp_filt.loc[f_filt, ratio_cols]

    ('plotting gaussian amplitude filter boxplots')
    fig, ax = plt.subplots(1, 1)
    nonsum_x = [0, 1.5, 3]
    gw_amp_x = [0.6, 2.1, 3.6]
    f_x = [1, 2.5, 4]

    IQR_boxplot(nonsum_id[ratio_cols], scatter_x=nonsum_x, color='k', ax=ax, label='NONSUM')
    IQR_boxplot(gw_amp_filt[ratio_cols], expecteds=expected_ratios, scatter_x=gw_amp_x, color='c', ax=ax,
                label='GW_AMP')
    IQR_boxplot(subf, scatter_x=f_x, color='m', ax=ax, label='AMP_F')
    plt.legend()
    fig.savefig(out + 'amps_filt_IQR.png')

    print('filtering by specified measure...')
    s_col = 'Gaussian_' + abbrev + '_s'
    w_col = 'Gaussian_' + abbrev + '_w'

    filt_peps = []
    for pep in gw_amp_filt.loc[f_filt, 'pep'].unique():
        pep_df = gw_amp_filt[gw_amp_filt['pep'] == pep]
        pep_s = fit_results.loc[pep_df.index[0], s_col]
        rts = pd.to_numeric(pep_df['retention_time'])
        rt_lb = rts.min()
        rt_ub = rts.max()

        if -1 * pep_s <= rt_ub and -1 * pep_s >= rt_lb:
            filt_peps.append(pep)

    sub_pepfilt = gw_amp_filt[gw_amp_filt['pep'].isin(filt_peps)]

    filt_scans = []
    for spectra in sub_pepfilt.index:
        rt = np.float(sub_pepfilt.loc[spectra, 'retention_time'])
        s = fit_results.loc[spectra, s_col]
        w = fit_results.loc[spectra, w_col]
        rt_lb = -1 * s - 0.75 * w
        rt_ub = -1 * s + 0.75 * w
        if rt >= rt_lb and rt <= rt_ub:
            filt_scans.append(spectra)

    sub_pepscanfilt = gw_amp_filt.loc[filt_scans]
    filtered_out = nonsum_id.loc[~nonsum_id.index.isin(filt_scans)]

    sub_pepscanfilt.to_csv(out + 'sub_pepscanfilt.csv')
    filtered_out.to_csv(out + 'filtered_out.csv')

    if args.plotfilt == True:
        os.mkdir(out + 'pepscanfilt/')
        os.mkdir(out + 'filtered_out/')

        print('plotting filtered and filtered-out spectra...')
        pepscanfilt_dfs = []
        for pep in sub_pepscanfilt['pep'].unique():
            sub = sub_pepscanfilt[sub_pepscanfilt['pep'] == pep]
            pepscanfilt_dfs.append([sub, out + 'pepscanfilt/', ft, wd])
        pool = Pool(args.threads)
        pool.map(plot_filt_spectra, pepscanfilt_dfs)

        filteredout_dfs = []
        for pep in filtered_out['pep'].unique():
            sub = filtered_out[filtered_out['pep'] == pep]
            filteredout_dfs.append([sub, out + 'filtered_out/', ft, wd])
        pool = Pool(args.threads)
        pool.map(plot_filt_spectra, filteredout_dfs)

        print('plotting filtered and filtered-out boxplots...')
        fig, ax = plt.subplots(1, 1)
        nonsum_x = [0, 1, 2]
        pepscanfilt_x = [0.25, 1.25, 2.25]
        filteredout_x = [0.5, 1.5, 2.5]
        IQR_boxplot(nonsum_id[ratio_cols], scatter_x=nonsum_x, color='k', ax=ax, label='Nonsum')
        IQR_boxplot(sub_pepscanfilt[ratio_cols], expecteds=expected_ratios, scatter_x=pepscanfilt_x, color='b', ax=ax,
                    label='Pepscanfilt')
        IQR_boxplot(filtered_out[ratio_cols], scatter_x=filteredout_x, color='g', ax=ax, label='Filtered_out')
        fig.savefig(out + 'pepscanfilt_IQR.png')

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    main(add_args(parser).parse_args())
