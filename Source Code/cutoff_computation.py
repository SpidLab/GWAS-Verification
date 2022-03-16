import csv
import pandas as pd
import argparse
import statsmodels.api as sm
import numpy as np
import math
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

parser = argparse.ArgumentParser(description="Parameters needed on the researcher side")
parser.add_argument('-i', '--dataset_file', type=str, help='Path to the original dataset file (E)', required=True)
parser.add_argument('-a', '--E_case_control_IDs_file', type=str, help='Path to the csv file that contains case and control IDs for dataset E', required=True)
parser.add_argument('-b', '--F_case_control_IDs_file', type=str, help='Path to the csv file that contains case and control IDs for dataset F', required=True)
parser.add_argument('-c', '--E_dataset_file', type=int, help='Path to the csv file that contains dataset E', required=True)
parser.add_argument('-d', '--F_dataset_file', type=int, help='Path to the csv file that contains dataset F', required=True)
parser.add_argument('-e', '--epsilon', type=int, help='Privacy parameter', required=True)
parser.add_argument('-s', '--shift', type=int, help='The index of SNP that is considered the most associated one. If s = 0, the correct SNPs are provided', required=True)
parser.add_argument('-l', '--l', type=int, help='Number of SNPs for which GWAS statistics are computed', required=True)
parser.add_argument('-n', '--nr_splits', type=int, help='Plot parameter to define number of split on the x-axis', required=True)
#If sufficient data is generated, increase nr_splits as it provides more accurate cutoff points
args = parser.parse_args()

F_case_control_IDs_file = "Data/F_case_control_IDs.csv"
E_case_control_IDs_file = "Data/E_case_control_IDs.csv"
F_dataset_file = "Data/F_original_dataset.csv"
E_dataset_file = "Data/E_original_dataset.csv"
epsilon = 3
shift = 0
l = 100
nr_splits = 10


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

def compute_relative_error(original_df, noisy_df, SNP_list):
    column_names = ['odds_ratio', 'MAF', 'p_val']
    RE_df = pd.DataFrame(columns=column_names)
    for SNP in SNP_list:
        odds_RE = MAF_RE = p_val_RE = 0.0000001
        if float(original_df.at[SNP, 'odds_ratio']) != 0:
            odds_RE = abs(float(original_df.at[SNP, 'odds_ratio']) - float(noisy_df.at[SNP, 'odds_ratio'])) / float(original_df.at[SNP, 'odds_ratio'])
        if float(original_df.at[SNP, 'MAF']) != 0:
            MAF_RE = abs(float(original_df.at[SNP, 'MAF']) - float(noisy_df.at[SNP, 'MAF'])) / float(original_df.at[SNP, 'MAF'])
        original_log_pval = -1*math.log(float(original_df.at[SNP, 'p_val']))
        noisy_log_pval = -1*math.log(float(noisy_df.at[SNP, 'p_val']))
        if original_log_pval != 0:
            p_val_RE = abs(original_log_pval-noisy_log_pval)/ original_log_pval
        RE_to_append = {"odds_ratio": odds_RE, 'MAF': MAF_RE, 'p_val': p_val_RE}
        RE_df = RE_df.append(RE_to_append, ignore_index=True)
    RE_df.index = SNP_list
    return RE_df

def compute_error(D_RE, E_RE, SNP_list):
    column_names = ['odds_ratio', 'MAF', 'p_val']
    error_df = pd.DataFrame(columns=column_names)
    for i in range(len(SNP_list)):
        odds_error = MAF_error = p_val_error = 0.0000001
        if E_RE.iat[i, 0] != 0:
            odds_error = abs(float(D_RE.iat[i, 0]) - float(E_RE.iat[i, 0])) /  float(E_RE.iat[i, 0])
        if E_RE.iat[i, 1] != 0:
            MAF_error = abs(float(D_RE.iat[i, 1]) - float(E_RE.iat[i, 1])) /  float(E_RE.iat[i, 1])
        if E_RE.iat[i, 2] != 0:
            p_val_error = abs(float(D_RE.iat[i, 2]) - float(E_RE.iat[i, 2])) /  float(E_RE.iat[i, 2])

        error_to_append = {"odds_ratio": odds_error, 'MAF': MAF_error, 'p_val': p_val_error}
        error_df = error_df.append(error_to_append, ignore_index=True)
    error_df.index = SNP_list
    return error_df

def generate_matlab_lines(pval_array):
    xaxis = []
    yaxis = []
    sum = 0
    max_val = math.ceil(max(pval_array))
    min_val = math.floor(min(pval_array))
    step = (max_val - min_val) / nr_splits
    for i in np.arange(min_val, max_val, step):
        xaxis.append(i)
        y = len([j for j in pval_array if j < i])
        diff = y - sum
        sum = y
        yaxis.append(diff / len(pval_array))

    file = open("Matlab/index.txt", "w")
    for i in xaxis:
        file.write(str(i) + " ")
    file.close()
    if (shift == 0):
        variable_name = "correct.txt"
    else:
        variable_name = "incorrect.txt"
    file = open("MATLAB/" + str(variable_name), "w")
    for i in yaxis:
        file.write(str(i) + " ")
    file.close()
    return 1

if __name__ == "__main__":
    # Get case user IDs and control user IDs of dataset F
    F_case_IDs, F_control_IDs = get_user_IDs(F_case_control_IDs_file)
    F_user_IDs = F_case_IDs + F_control_IDs

    # Get case user IDs and control user IDs of dataset E
    E_case_IDs, E_control_IDs = get_user_IDs(E_case_control_IDs_file)
    E_user_IDs = E_case_IDs + E_control_IDs

    # Load dataset F
    F_df = pd.read_csv(F_dataset_file, sep=',', index_col=0)
    print("Dataset F loaded successfully.")

    # Perform GWAS on dataset F
    F_GWAS = perform_GWAS(F_df, F_case_IDs, F_control_IDs, False)
    print("GWAS performed on dataset F.")

    # Sort according to p-value
    F_GWAS = F_GWAS.sort_values(by='p_val', ascending=True)  # sort according to p value in ascending order
    original_SNPs = F_GWAS.index[0:l]  # keep track of the original SNP IDs, we only shift the values
    F_GWAS = F_GWAS[shift:shift + l]  # pick the top l SNPs
    F_GWAS.index = original_SNPs  # reassign the original IDs, but with other statistics
    print("GWAS performed on dataset E.")
    SNP_list = F_GWAS.index

    # Truncate F by selecting only the SNPs that are provided in GWAS
    F_df = F_df[F_df.index.isin(SNP_list)]
    F_df = F_df.reindex(SNP_list)

    # Add noise to F
    F_noisy_df = generate_noisy_dataframe(F_df, F_user_IDs)

    # Perform GWAS on noisy dataset F without using aggreagation
    F_noisy_GWAS = perform_GWAS(F_noisy_df, F_case_IDs, F_control_IDs, False)


    # Load dataset E
    E_df = pd.read_csv(E_dataset_file, sep=',', index_col=0)
    print("Dataset E loaded successfully.")

    # Perform GWAS on the original dataset E
    E_GWAS = perform_GWAS(E_df, E_case_IDs, E_control_IDs, False)
    sorted_E_GWAS = E_GWAS.sort_values(by='p_val', ascending=True)  # sort according to p value in ascending order
    E_GWAS = sorted_E_GWAS[0:l]
    print("GWAS performed on dataset E.")

    # Add noise to E
    E_df = E_df.reindex(E_GWAS.index)  # reindex the dataframe according to p value
    E_noisy_df = generate_noisy_dataframe(E_df[0:0 + l], E_user_IDs)

    # Perform GWAS on noisy dataset E without using aggreagation
    E_noisy_GWAS = perform_GWAS(E_noisy_df, E_case_IDs, E_control_IDs, False)

    # Compute the deviation - relative error(RE) for each of the statistics of dataset D, E and store in dataframes D_RE and E_RE respectively
    F_RE = compute_relative_error(F_GWAS, F_noisy_GWAS, SNP_list)
    E_SNP_list = E_GWAS.index
    E_RE = compute_relative_error(E_GWAS, E_noisy_GWAS, E_SNP_list)

    # Compute the error between deviations
    error_df = compute_error(F_RE, E_RE, SNP_list)

    pval_array = error_df['p_val'].tolist()
    generate_matlab_lines(pval_array)
    print(pval_array)

    #Similarly, we can compute the cutoff point for odds ratio and MAF