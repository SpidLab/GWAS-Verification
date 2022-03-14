import csv
import pandas as pd
import argparse
import statsmodels.api as sm
import numpy as np
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

# parser = argparse.ArgumentParser(description="Parameters needed on the researcher side")
# parser.add_argument('-i', '--dataset_dir', type=str, help='Path to the original dataset file', required=True)
# parser.add_argument('-e', '--epsilon', type=int, help='Privacy parameter.', required=True)
# parser.add_argument('-k', '--k', type=int, help='Number of SNPs that are provided in the partial noisy dataset', required=True)
# parser.add_argument('-l', '--l', type=int, help='Number of SNPs for which GWAS statistics are provided', required=True)
# parser.add_argument('-s', '--shift', type=int, help='The index of SNP that is considered the most associated one. If s = 0, the correct SNPs are provided', required=True)
# parser.add_argument('-c', '--case_control_IDs_file', type=str, help='Path to the csv files that contains case and control IDs', required=True)
# parser.add_argument('-o', '--output_dir', type=str, help='Path to save the output files.', required=True)
# args = parser.parse_args()

#parameters
dataset_file = "Data/E_original_dataset.csv"
E_case_control_IDs_file = "Data/E_case_control_IDs.csv"
noisy_dataset_file = "Data/Researcher/noisy_dataset_D.csv"
D_dataset_GWAS_file = "Data/Researcher/D_GWAS.csv"
D_case_control_IDs_file = "Data/D_case_control_IDs.csv"
epsilon = 3
odds_threshold = 1.1
maf_threshold = 1.2
pval_threshold = 1.3
#todo change thresholds to a file
#output_dir = "Data/Verifier/"

def get_user_IDs(case_control_IDs_file):
    with open(case_control_IDs_file) as csvfile:
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

def estimate_value(dataframe, state):
    p = np.exp(epsilon) / (np.exp(epsilon) + 2)
    q = 1 / (np.exp(epsilon) + 2)
    n = len(dataframe.columns)
    if (p != q):
        ci = dataframe.apply(lambda x: (x == state).sum(), axis=1)
        cpi = ((ci - n * q) / (p - q))
        cpi = cpi.apply(np.int64)
        return cpi
    else:
        return dataframe.apply(lambda x: (x == state).sum(), axis=1)


#if aggregation technique is used, the verifier estimates the actual occurrences of 0, 1, and 2
def perform_GWAS(dataframe, case_IDs, control_IDs, aggregation):
    case_dataframe = dataframe[case_IDs]
    control_dataframe = dataframe[control_IDs]

    if aggregation == False:
        dataframe['case_0'] = case_dataframe.apply(lambda x: (x == 0).sum(), axis=1)
        dataframe['case_1'] = case_dataframe.apply(lambda x: (x == 1).sum(), axis=1)
        dataframe['case_2'] = case_dataframe.apply(lambda x: (x == 2).sum(), axis=1)
        dataframe['control_0'] = control_dataframe.apply(lambda x: (x == 0).sum(), axis=1)
        dataframe['control_1'] = control_dataframe.apply(lambda x: (x == 1).sum(), axis=1)
        dataframe['control_2'] = control_dataframe.apply(lambda x: (x == 2).sum(), axis=1)
    else:
        dataframe['case_0'] = estimate_value(case_dataframe, 0)
        dataframe['case_1'] = estimate_value(case_dataframe, 1)
        dataframe['case_2'] = estimate_value(case_dataframe, 2)
        dataframe['control_0'] = estimate_value(control_dataframe, 0)
        dataframe['control_1'] = estimate_value(control_dataframe, 1)
        dataframe['control_2'] = estimate_value(control_dataframe, 2)

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
    p = np.exp(epsilon) / (np.exp(epsilon) + 2)
    q = 1 / (np.exp(epsilon) + 2)
    for j in range(len(user_IDs)):
        dataframe[str(user_IDs[j])] = dataframe.apply(lambda x: randomized_response(x[str(user_IDs[j])], p, q), axis=1)
    return  dataframe


if __name__ == "__main__":
    # Get case user IDs and control user IDs of dataset D
    D_case_IDs, D_control_IDs = get_user_IDs(D_case_control_IDs_file)
    D_user_IDs = D_case_IDs + D_control_IDs

    # Get case user IDs and control user IDs of dataset E
    E_case_IDs, E_control_IDs = get_user_IDs(E_case_control_IDs_file)
    E_user_IDs = E_case_IDs + E_control_IDs

    # Load GWAS of dataset D
    D_GWAS = pd.read_csv(D_dataset_GWAS_file, sep=',', index_col=0)
    print("GWAS of dataset D loaded successfully.")

    # Load partial noisy dataset D
    D_noisy_df = pd.read_csv(noisy_dataset_file, sep =',', index_col=0)
    D_noisy_df = D_noisy_df[D_noisy_df.index.isin(D_GWAS.index)] # get only the SNPs that are provided in GWAS stats
    D_noisy_df = D_noisy_df.reindex(D_GWAS.index)
    print("Noisy dataset D loaded successfully.")

    # Perform GWAS on noisy dataset D using aggreagation
    D_noisy_GWAS = perform_GWAS(D_noisy_df, D_case_IDs, D_control_IDs, True)
    print("GWAS performed on noisy dataset D.")

    # Load dataset E
    E_df = pd.read_csv(dataset_file, sep=',', index_col=0)
    print("Dataset E loaded successfully.")

    # Truncate E by selecting only the SNPs that are provided in GWAS of D
    E_df = E_df[E_df.index.isin(D_GWAS.index)]
    E_df = E_df.reindex(D_GWAS.index)

    # Perform GWAS on the original dataset E
    E_GWAS = perform_GWAS(E_df, E_case_IDs, E_control_IDs, False)
    print("GWAS performed on dataset E.")

    # Add noise to E
    E_noisy_df = generate_noisy_dataframe(E_df, E_user_IDs)

    # Perform GWAS on noisy dataset E using aggreagation
    E_noisy_GWAS = perform_GWAS(E_noisy_df, E_case_IDs, E_control_IDs, True)

    # Include the code that computes RE
    # maf_error = []
    # mre_or_error = []
    # pval_error = []
    # for k in range(end_snp):
    #     maf_error.append(abs(or_stats.aaf.iloc[k] - noisy_mix_res.aaf.iloc[k]) / or_stats.aaf.iloc[k])
    #     mre_or_error.append(
    #         abs(or_stats.odds_ratio.iloc[k] - noisy_mix_res.odds_ratio.iloc[k]) / or_stats.odds_ratio.iloc[k])
    #     pval_error.append(abs(or_stats.log_pval.iloc[k] - noisy_mix_res.log_pval.iloc[k]) / or_stats.log_pval.iloc[k])