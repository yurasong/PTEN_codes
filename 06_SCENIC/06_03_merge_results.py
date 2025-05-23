#######################################################################################
# 06_03_merge_results.py
# Please activate conda environment before running this file.
# Please run this on python environment.
# This code will return average value of AUC score per regulon as a matrix.
# Acknowledgement: LEE Jeong Gyu (Hanyang Univ. College of Communications, Information Sociology, leejg794294@gmail.com)
#######################################################################################


import pandas as pd
from glob import glob

# 1) Gather all auc_mtx files
file_paths = sorted(glob("auc_mtx*.csv"))  # matches auc_mtx.csv, auc_mtx1.csv, …, auc_mtx9.csv

# 2) Read each into a DataFrame, keeping Cell as a column
dfs = [pd.read_csv(fp) for fp in file_paths]

# 3) Use the first DF's Cell column as our index
result = pd.DataFrame(dfs[0]["Cell"], columns=["Cell"])

# 4) Concatenate all on Cell, creating duplicate column names where they overlap
concat = pd.concat(
    [df.set_index("Cell") for df in dfs],
    axis=1,
    sort=False
)

# 5) Group-by column name and take the mean across the duplicates
averaged = (
    concat
    .groupby(level=0, axis=1)  # groups columns by their name
    .mean()
    .reset_index()
)

# 6) Write out
averaged.to_csv("AUC_Average.csv", index=False)
