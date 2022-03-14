import csv
import pandas as pd
import argparse
import statsmodels.api as sm
import numpy as np
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

parser = argparse.ArgumentParser(description="Generate dataset")
parser.add_argument('-i', '--dataset_dir', type=str, help='Path to the original dataset file', required=True)
parser.add_argument('-e', '--epsilon', type=int, help='Privacy parameter.', required=True)
parser.add_argument('-k', '--k', type=int, help='Number of SNPs that are provided in the partial noisy dataset', required=True)
parser.add_argument('-l', '--l', type=int, help='Number of SNPs for which GWAS statistics are provided', required=True)
parser.add_argument('-s', '--shift', type=int, help='The index of SNP that is considered the most associated one. If s = 0, the correct SNPs are provided', required=True)
parser.add_argument('-c', '--case_control_IDs_file', type=str, help='Path to the csv files that contains case and control IDs', required=True)
parser.add_argument('-o', '--output_dir', type=str, help='Path to save the output files.', required=True)
args = parser.parse_args()

def get_user_IDs():
    with open(args.case_control_IDs_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for index, row in enumerate(spamreader):
            if index == 0:
                case_IDs = row
            if index == 1:
                control_IDs = row
    return case_IDs, control_IDs

def compute_statistics(a, b, c, d):
    table = sm.stats.Table2x2(np.array([[a, b], [c, d]]))
    low_interval, high_interval = table.oddsratio_confint(alpha=0.05, method='normal')
    column_names = ['odds_ratio', 'low_interval', 'high_interval', 'se', 'p_val']
    return pd.Series(["%.8f" % table.oddsratio, "%.8f" % low_interval, "%.8f" % high_interval, "%.8f" % table.log_oddsratio_se, "%.8f" % table.log_oddsratio_pvalue()],  index=column_names)


def perform_GWAS(dataframe, case_IDs, control_IDs):
    case_dataframe = dataframe[case_IDs]
    control_dataframe = dataframe[control_IDs]

    dataframe['case_0'] = case_dataframe.apply(lambda x: (x == 0).sum(), axis=1)
    dataframe['case_1'] = case_dataframe.apply(lambda x: (x == 1).sum(), axis=1)
    dataframe['case_2'] = case_dataframe.apply(lambda x: (x == 2).sum(), axis=1)
    dataframe['control_0'] = control_dataframe.apply(lambda x: (x == 0).sum(), axis=1)
    dataframe['control_1'] = control_dataframe.apply(lambda x: (x == 1).sum(), axis=1)
    dataframe['control_2'] = control_dataframe.apply(lambda x: (x == 2).sum(), axis=1)

    dataframe['case_major'] = dataframe['case_0']
    dataframe['case_minor'] = dataframe['case_1'] + dataframe['case_2']
    dataframe['case_minor_counts'] = dataframe['case_1'] + 2 * dataframe['case_2']
    dataframe['control_major'] = dataframe['control_0']
    dataframe['control_minor'] = dataframe['control_1'] + dataframe['control_2']
    dataframe['control_minor_counts'] = dataframe['control_1'] + 2 * dataframe['control_2']

    dataframe_GWAS = dataframe.apply(
    lambda x: compute_statistics(x['case_minor'], x['control_minor'], x['case_major'], x['control_major']), axis=1)
    dataframe_GWAS['alt'] = dataframe['alt']
    dataframe_GWAS['ref'] = dataframe['ref']
    dataframe_GWAS['case_MAF'] = dataframe['case_minor_counts'] / (2.0 * len(case_IDs))
    dataframe_GWAS['control_MAF'] = dataframe['control_minor_counts'] / (2.0 * len(control_IDs))
    dataframe_GWAS['MAF'] = (dataframe['case_minor_counts'] + dataframe['control_minor_counts']) / (
                2.0 * (len(case_IDs) + len(control_IDs)))
    return dataframe_GWAS

def randomized_response(val, p, q):
    rand_val = np.random.uniform(0, 1)
    new_val = val
    if rand_val > p:
        if rand_val > p + q:
            if val == 1:
                new_val = 2
            else:
                new_val = 2 - val
        else:
            new_val = abs(1 - val)
    return new_val

def generate_noisy_dataframe(dataframe, user_IDs):
    p = np.exp(args.epsilon) / (np.exp(args.epsilon) + 2)
    q = 1 / (np.exp(args.epsilon) + 2)
    for j in range(len(user_IDs)):
        # dataframe.loc[:,str(user_IDs[j])] = dataframe.apply(lambda x: randomized_response(x[str(user_IDs[j])], p, q), axis=1)
        dataframe[str(user_IDs[j])] = dataframe.apply(lambda x: randomized_response(x[str(user_IDs[j])], p, q), axis=1)
    return  dataframe


if __name__ == "__main__":
    # Get case user IDs and control user IDs
    case_IDs, control_IDs = get_user_IDs()
    user_IDs = case_IDs + control_IDs
    print("The total number of samples is: " + str(len(user_IDs)))

    # Load dataset D
    D_df = pd.read_csv(args.dataset_dir, sep =',', index_col=0)
    D_df = D_df.head(200) #todo remove this line - kept it 200 for testing
    print("Dataset D loaded successfully.")

    # Perform GWAS on the original dataset D
    D_GWAS = perform_GWAS(D_df, case_IDs, control_IDs)
    print("GWAS performed on dataset D.")

    sorted_D_GWAS = D_GWAS.sort_values(by='p_val', ascending=True) #sort according to p value in ascending order
    shared_D_GWAS = sorted_D_GWAS[args.shift:args.shift + args.l] #pick the top l SNPs
    shared_D_GWAS.to_csv(args.output_dir + "D_GWAS.csv")
    print("Generated GWAS results for top " + str(args.l) + " SNPs.")

    # pick top k SNPs from sorted_D_GWAS
    sorted_D_df = D_df.reindex(sorted_D_GWAS.index) #reindex the dataframe according to p value
    print(sorted_D_df)
    noisy_D = generate_noisy_dataframe(sorted_D_df[args.shift:args.shift+args.k], user_IDs) # add noise to only k SNPs
    noisy_D.to_csv(args.output_dir + "noisy_dataset_D.csv")
    print("Generated partial noisy dataset for top " + str(args.k) + " SNPs and " + str(len(user_IDs)) + " users.")
