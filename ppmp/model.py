#!/usr/bin/env python3

import glob
import json
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns
from pandas.plotting import parallel_coordinates
from scipy import stats
from tqdm import tqdm

from ppmp.protein import Protein
from ppmp.scripts.sci_utils import is_outlier
from ppmp.scripts.utils import alpha_label
from ppmp.scripts.utils import handle_missing_dirs


def analysis(csv_path='./data/csv/',
             test_path=None, sim=True, parallel=True, dist=True, violin=True, whitney=True, order=True, lvplot=True):
    protein_df = read_csv(csv_path=csv_path)
    all_df = pd.concat(protein_df)

    if parallel:
        parallel_plot(protein_df=protein_df)

    all_df = create_triplets(all_df=all_df)

    if dist:
        dist_plot(all_df=all_df)

    if violin:
        violin_plot(all_df)

    if lvplot:
        lv_plot(all_df)

    triplet_df = mann_whitney(all_df, whitney)
    w = list()
    p = list()
    normal = list()
    mean = list()
    std = list()
    for triplet in tqdm(triplet_df['triplet'],
                        desc='Shapiro normal test statistics and mean/standard deviation for triplet'):
        triplet_rmsd = all_df[all_df['triplet'] == triplet]['rmsd']
        filtered = triplet_rmsd[~is_outlier(triplet_rmsd)]
        shapiro_result = shapiro_test(filtered)
        w.append(shapiro_result[0])
        p.append(shapiro_result[1])
        normal.append(shapiro_result[2])
        mean.append(np.mean(filtered))
        std.append(np.std(filtered))
    triplet_df['w'] = w
    triplet_df['shapiro p'] = p
    triplet_df['normal'] = normal
    triplet_df['mean'] = mean
    triplet_df['std'] = std

    single_df = pd.DataFrame(columns=['module', 'mean', 'std'])
    for single in tqdm(all_df['module'].unique(), desc='Mean/standard deviation for single'):
        rmsd_single = all_df[all_df['module'] == single]['rmsd']
        single_df.loc[single_df.shape[0]] = [single,
                                             np.mean(rmsd_single),
                                             np.std(rmsd_single)]
    # Validate that order matters
    if order:
        ordered_df = _validate_order_significance(all_df)
        print(ordered_df[ordered_df['p'] < 0.05].shape)

    if test_path is not None:
        test_protein = dict()
        if os.path.basename(test_path) == '':
            for path in tqdm(glob.glob(test_path + '*.json'), desc='Reading the test JSON into the memory'):
                file_name = os.path.basename(path)
                protein_name = os.path.splitext(file_name)[0]
                with open(path, 'r') as f:
                    test_protein[protein_name] = json.load(f)

        else:
            print('Reading the test JSON into the memory')
            file_name = os.path.basename(test_path)
            protein_name = os.path.splitext(file_name)[0]
            with open(test_path, 'r') as f:
                test_protein[protein_name] = json.load(f)

        all_test = dict()
        for protein_name, protein_data in tqdm(test_protein.items(), desc='Creating triplets for test data'):
            test_df = pd.DataFrame()
            test_df['module'] = protein_data['nodes']
            test_triplet = list()
            for i, v in enumerate(protein_data['nodes']):
                test_triplet.append(
                    (
                        _access_list(i - 1, protein_data['nodes']),
                        v,
                        _access_list(i + 1, protein_data['nodes'])
                    )
                )
            test_df['triplet'] = test_triplet
            all_test[protein_name] = test_df

        if sim:
            for protein_name in tqdm(all_test, desc='Calculating RMSD for the test data'):

                original = Protein(protein_name,
                                   test_path + protein_name + '.pdb',
                                   test_path + protein_name + '.json',
                                   strict=True)

                temp_rmsd = dict()
                for path in glob.glob('{}{}_*.pdb'.format(test_path, protein_name)):
                    variation = os.path.splitext(os.path.basename(path))[0]

                    perturbed = Protein(variation,
                                        path,
                                        test_path + protein_name + '.json',
                                        strict=True)
                    temp_rmsd[variation] = Protein.kabsch(original, perturbed, range(len(original.modules_chain)))

                temp_pro = list()
                for order, mod in enumerate(all_test[protein_name]['module']):
                    temp_var = np.mean([temp_rmsd[var][order] for var in temp_rmsd])
                    temp_pro.append(temp_var)

                all_test[protein_name]['rmsd mean'] = temp_pro

        # Make a prediction
        for protein_name in tqdm(all_test, desc='Making and exporting predictions'):
            all_test[protein_name]['single dist'] = all_test[protein_name]['module'].apply(single_dist,
                                                                                           args=(single_df,))
            all_test[protein_name]['prediction'] = all_test[protein_name]['triplet'].apply(predict_module,
                                                                                           args=(single_df, triplet_df))
            all_test[protein_name]['triplet pred'] = \
                all_test[protein_name]['single dist'] != all_test[protein_name]['prediction']
            all_test[protein_name] = _prediction_to_csv(all_test[protein_name], protein_name)
            all_test[protein_name]['abs pred err'] = \
                abs(all_test[protein_name]['rmsd mean'] - all_test[protein_name]['prediction mean'])

        for protein_name in tqdm(all_test, desc='Exporting and plotting predictions'):
            x = np.linspace(0, 5, 10000)
            for i, r in all_test[protein_name].iterrows():
                single_normal = stats.norm.pdf(x, r['single mean'], r['single std'])
                plt.plot(x, single_normal, label=r['module'])
                if r['triplet pred']:
                    triplet_normal = stats.norm.pdf(x, r['prediction mean'], r['prediction std'])
                    plt.plot(x, triplet_normal, label=r['triplet'])

                plt.axvline(r['rmsd mean'], color='r', label='Actual mean of RMSD')
                plt.xlabel('RMSD')
                plt.ylabel('Frequency')
                plt.legend()
                figs_path = os.path.join(os.getcwd(), 'out', 'figs', 'prediction', protein_name + '-' + str(i) + '.pdf')
                plt.savefig(handle_missing_dirs(figs_path))
                plt.close()

        test_table = pd.concat([i for i in all_test.values()])
        correct = test_table['prediction mean'] - test_table['prediction std'] <= test_table['rmsd mean']
        correct = correct & (test_table['rmsd mean'] <= test_table['prediction mean'] + test_table['prediction std'])
        test_table = test_table[~ correct].reset_index(drop=True)

        x = np.arange(len(test_table))
        myx_ticks = test_table['module']

        plt.figure(figsize=(17, 14))
        plt.xticks(x, myx_ticks)
        plt.title('Plot Demonstrating the Error of the Incorrectly Predicted Modules', fontsize=20)
        plt.xlabel('Module', fontsize=20)
        plt.ylabel('RMSD', fontsize=20)

        for i in tqdm(range(len(test_table['module'])), desc='Plotting error graph'):
            if test_table['triplet pred'][i]:
                plt.scatter(x[i], test_table['prediction mean'][i], c='b', marker='o', s=75)
                plt.errorbar(x[i],
                             test_table['prediction mean'][i],
                             yerr=test_table['prediction std'][i],
                             fmt='o',
                             c='b',
                             elinewidth=2)
            else:
                plt.scatter(x[i], test_table['prediction mean'][i], c='r', marker='o', s=75)
                plt.errorbar(x[i],
                             test_table['prediction mean'][i],
                             yerr=test_table['prediction std'][i],
                             fmt='o',
                             c='r',
                             elinewidth=2)

        plot = plt.scatter(x, test_table['rmsd mean'], marker='o', c='black', label='The average RMSD from simulation')
        red_patch = mpatches.Patch(color='r', label='Predictions made from a module distribution')
        blue_patch = mpatches.Patch(color='b', label='Predictions made from the triplet distribution')
        data_patch = mpatches.Patch(color='black', label='The average RMSD from simulation')
        plt.legend(handles=[red_patch, blue_patch, plot], prop={'size': 20})

        plt.savefig("./out/prediction/analysis.pdf")
        plt.close()

        return all_df, single_df, triplet_df, all_test
    else:
        return all_df, single_df, triplet_df


def read_csv(csv_path):
    """Read CSV into dictionary of DataFrames

    :param csv_path: Path to the folder containing CSV files to be read.`
    :return: Dictionary of DataFrames. One CSV file translates to one DataFrame.
    """
    protein_df = dict()

    for path in tqdm(glob.glob(csv_path + '*.csv'),
                     desc='Reading the CSV into the memory'):
        file_name = os.path.basename(path)
        protein_name = os.path.splitext(file_name)[0]
        protein_df[protein_name] = pd.read_csv(path)

    return protein_df


def parallel_plot(protein_df):
    """Parallel plot per protein.

    :param protein_df: Dictionary of proteins/DataFrames. One protein is represented by one DataFrame.
    :return: Save the parallel plots to the disk.
    """
    for protein_name in tqdm(protein_df, desc='Creating parallel plots'):
        pivot_df = protein_df[protein_name].pivot_table('rmsd', ['order', 'module'], 'variation')
        pivot_df['module'] = pivot_df.index.get_level_values('module')

        parallel_coordinates(pivot_df, 'module', colormap=plt.get_cmap("Set2"), linewidth=5)  # TODO Better output
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.xticks(rotation=60)
        plt.xlabel('Perturbations')
        plt.ylabel('RMSD')

        plt.savefig(handle_missing_dirs('./out/figs/parallel/' + protein_name + '.pdf'), bbox_inches="tight")
        # plt.savefig('./out/figs/parallel/png/' + protein_name + '.png', bbox_inches="tight")
        plt.close()


def create_triplets(all_df):
    """Create triplet.

    :param all_df: DataFrame for all proteins/modules.
    :return: DataFrame with triplets included.
    """
    triplet_list = list()
    # unordered_list = list()
    for variation in tqdm(all_df['variation'].unique(), desc='Creating triplets for training data'):

        module_list = list(all_df[all_df['variation'] == variation]['module'])
        order_list = list(all_df[all_df['variation'] == variation]['order'])

        for i, v in zip(order_list, module_list):
            triplet_list.append(
                (
                    _access_list(i - 1, module_list),
                    _access_list(i, module_list),
                    _access_list(i + 1, module_list)
                )
            )

    all_df['triplet'] = triplet_list
    # all_df['unordered'] = unordered_list
    return all_df


def dist_plot(all_df):

    for single in tqdm(all_df['module'].unique(), desc='Creating histogram/KDE plots'):
        f, (ax_box, ax_hist) = plt.subplots(2, sharex='all', gridspec_kw={"height_ratios": (.15, .85)})

        single_rmsd_array = all_df[all_df['module'] == single]['rmsd']

        temp_df = all_df[all_df['module'] == single]
        for triplet in temp_df['triplet'].unique():
            triplet_rmsd_list = temp_df[temp_df['triplet'] == triplet]['rmsd']
            sns.kdeplot(triplet_rmsd_list, label=triplet, ax=ax_hist)

        ax_box.set(xlabel='')
        ax_hist.legend().set_visible(False)
        sns.boxplot(single_rmsd_array, ax=ax_box)
        sns.distplot(single_rmsd_array, ax=ax_hist, kde=True, hist=True, rug=True)
        f.suptitle(single)
        f.savefig(handle_missing_dirs('./out/figs/kde-hist/' + single + '.pdf'))
        plt.close(f)


def violin_plot(all_df):
    for single in tqdm(all_df['module'].unique(), desc='Creating violin plots'):
        plt.xticks(rotation=90)
        xlabel, lookup = alpha_label(all_df[all_df['module'] == single]['triplet'])
        sns.violinplot(x=xlabel,
                       y=all_df[all_df['module'] == single]['rmsd'])
        plt.tight_layout()
        plt.savefig(handle_missing_dirs('./out/figs/violin/' + single + '.pdf'))  # TODO Better output
        plt.close()


def lv_plot(all_df):
    sns.lvplot(x='module', y='rmsd', data=all_df.sort_values(['rmsd'], ascending=[True]), scale='linear', palette='mako')
    plt.xticks(rotation=60)
    plt.suptitle('RMSD Range per Module')
    plt.xlabel('Module')
    plt.ylabel('RMSD')
    plt.tight_layout()
    plt.savefig(handle_missing_dirs('./out/figs/lv_plot.pdf'))


def _access_list(index, iterable):
    """Accessing list by index, different behaviour for negative or
    out-of-bounds indices.

    :param index: An index we wan to access.
    :param iterable: An indexable iterable object which we want to access.
    :return: Return the stored value for the corresponding index, if the
             index is negative or out-of-bounds, return 'EMPTY' string.
    """
    if index < 0:
        return 'EMPTY'
    else:
        try:
            return iterable[index]
        except IndexError:
            return 'EMPTY'


def mann_whitney(all_df, plot):
    from statsmodels.stats import multitest
    df = pd.DataFrame(columns=['triplet', 'u', 'p'])
    for triplet in tqdm(all_df['triplet'].unique(), desc='Calculating Mann-Whitney test statistic'):
        x = all_df[all_df['module'] == triplet[1]]['rmsd'].values
        y = all_df[all_df['triplet'] == triplet]['rmsd'].values
        u, p = stats.mannwhitneyu(x, y)
        if plot:
            sns.distplot(x, rug=True)
            sns.distplot(y, rug=True)
            plt.suptitle('{}\nU statistic: {}    p-value: {}'.format(triplet, u, p))
            plt.savefig(handle_missing_dirs('./out/figs/mann_whitney/' + str(triplet) + '.pdf'))
            plt.close()
        df.loc[df.shape[0]] = [triplet, u, p]

    correction = multitest.multipletests(df['p'], alpha=0.05, method='hs')
    df['reject'] = correction[0]
    df['p-value corrected'] = correction[1]
    return df


def shapiro_test(rmsd):
    w, p = stats.shapiro(rmsd)
    normal = p > 0.05
    return w, p, normal


def predict_module(triplet, single_df, triplet_df):

    triplet_row = triplet_df[triplet_df['triplet'] == triplet]

    predicted_mean = single_df[single_df['module'] == triplet[1]]['mean'].values[0]
    predicted_std = single_df[single_df['module'] == triplet[1]]['std'].values[0]

    if not triplet_row.empty:
        try:
            reject = triplet_row['reject'].all()
            normal = triplet_row['normal'].all()
        except ValueError:
            reject = False
            normal = False

        if reject and normal:
            predicted_mean = triplet_row['mean'].values[0]
            predicted_std = triplet_row['std'].values[0]

    return predicted_mean, predicted_std


def single_dist(single, single_df):

    mean = single_df[single_df['module'] == single]['mean'].values[0]
    std = single_df[single_df['module'] == single]['std'].values[0]

    return mean, std


def _prediction_to_csv(protein_df, protein_name):
    new_col_list = ['prediction mean', 'prediction std']
    for n, col in enumerate(new_col_list):
        protein_df[col] = protein_df['prediction'].apply(lambda prediction: prediction[n])

    protein_df = protein_df.drop('prediction', axis=1)

    new_col_list = ['single mean', 'single std']
    for n, col in enumerate(new_col_list):
        protein_df[col] = protein_df['single dist'].apply(lambda prediction: prediction[n])

    protein_df = protein_df.drop('single dist', axis=1)

    directory = os.path.join(os.getcwd(), 'out', 'prediction')
    protein_df.to_csv(handle_missing_dirs(os.path.join(directory, protein_name + '.csv')))

    return protein_df


def _validate_order_significance(all_df):
    from statsmodels.stats import multitest
    unordered_df = pd.DataFrame(columns=['A', 'B', 'u', 'p'])
    for A in all_df['triplet'].unique():
        B = A[::-1]
        if A not in unordered_df['A'].unique() and A not in unordered_df['B'].unique():
            rmsd_list_A = all_df[all_df['triplet'] == A]['rmsd']
            rmsd_list_B = all_df[all_df['triplet'] == B]['rmsd']
            if not (rmsd_list_A.empty or rmsd_list_B.empty):
                u, p = stats.mannwhitneyu(rmsd_list_A, rmsd_list_B)
                unordered_df.loc[unordered_df.shape[0]] = [A, B, u, p]
    correction = multitest.multipletests(unordered_df['p'], alpha=0.05, method='hs')
    unordered_df['reject'] = correction[0]
    unordered_df['p-value corrected'] = correction[1]

    return unordered_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('csv_path', type=str, help='The folder containing the CSV files for training.')
    parser.add_argument('-t', '--test', type=str, help='The folder (or file) containing data we want to predict.')
    parser.add_argument('-p', '--parallel', action='store_true', help='Create parallel plot.')
    parser.add_argument('-v', '--violin', action='store_true', help='Create violin plot.')
    parser.add_argument('-d', '--dist', action='store_true', help='Create distribution plot.')
    parser.add_argument('-w', '--whitney', action='store_true', help='Whitney test')
    parser.add_argument('-s', '--sim', action='store_true', help='Compare prediction with simulation')
    parser.add_argument('-o', '--order', action='store_true', help='Validate assumption that order matters.')
    parser.add_argument('-l', '--lvplot', action='store_true', help='Create lvplot.')
    args = parser.parse_args()
    analysis(csv_path=args.csv_path,
             test_path=args.test_path,
             parallel=args.parallel,
             violin=args.violin,
             dist=args.dist,
             whitney=args.whitney,
             sim=args.sim,
             order=args.order,
             lvplot=args.lvplot)
