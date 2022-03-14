import csv
import pandas as pd
import argparse
import statsmodels.api as sm
import numpy as np

#List of input parameters - previous method
dataset_path = "Data/original_dataset_D.csv"
epsilon = 3
k = 150 # compute noisy D based on K
l = 100 # provide GWAS for l SNPs
flag = "Correct"
user_IDs_file = "Data/case_control_IDs.csv"

# parser = argparse.ArgumentParser(description="Generate dataset")
# parser.add_argument('-i', '--data_path', type=str, help='Path to the directory which contains the users (samples).', required=True)
# parser.add_argument('-c', '--chromosome_number', type=str, help='Pick for which chromosome will be extracted the data.', required=True)
# parser.add_argument('-s', '--sequenced_user_ID', type=str, help='Path to the user from which the SNP IDs will be selected.', required=True)
# parser.add_argument('-u', '--user_IDs_file', type=str, help='Path to the csv files that contains case and control IDs.', required=True)
# parser.add_argument('-o', '--output_file', type=str, help='Path to save the output file.', required=True)
# args = parser.parse_args()

def get_user_IDs():
    with open(user_IDs_file) as csvfile:
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
    dataframe_GWAS = dataframe.apply(
    lambda x: compute_statistics(x['case_minor'], x['control_minor'], x['case_major'], x['control_major']), axis=1)
    dataframe_GWAS['alt'] = dataframe['alt']
    dataframe_GWAS['ref'] = dataframe['ref']
    dataframe_GWAS['case_aaf'] = dataframe['case_minor_counts'] / (2.0 * len(case_IDs))
    dataframe_GWAS['control_aaf'] = dataframe['control_minor_counts'] / (2.0 * len(control_IDs))
    dataframe_GWAS['aaf'] = (dataframe['case_minor_counts'] + dataframe['control_minor_counts']) / (
                2.0 * (len(case_IDs) + len(control_IDs)))
    return dataframe_GWAS

if __name__ == "__main__":
    # Get case user IDs and control user IDs
    case_IDs, control_IDs = get_user_IDs()
    user_IDs = case_IDs + control_IDs
    print("The total number of samples is: " + str(len(user_IDs)))

    D_df = pd.read_csv(dataset_path, sep =',', index_col=0)
    print(D_df)

    D_GWAS = perform_GWAS(D_df, case_IDs, control_IDs)
    print(D_GWAS)