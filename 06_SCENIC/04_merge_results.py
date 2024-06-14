import pandas as pd

## Import our data by panda

mtx = pd.read_csv("auc_mtx.csv")
mtx1 = pd.read_csv("auc_mtx1.csv")
mtx2 = pd.read_csv("auc_mtx2.csv")
mtx3 = pd.read_csv("auc_mtx3.csv")
mtx4 = pd.read_csv("auc_mtx4.csv")
mtx5 = pd.read_csv("auc_mtx5.csv")
mtx6 = pd.read_csv("auc_mtx6.csv")
mtx7 = pd.read_csv("auc_mtx7.csv")
mtx8 = pd.read_csv("auc_mtx8.csv")
mtx9 = pd.read_csv("auc_mtx9.csv")

mtx_columns = list(mtx.columns)
mtx1_columns = list(mtx1.columns)
mtx2_columns = list(mtx2.columns)
mtx3_columns = list(mtx3.columns)
mtx4_columns = list(mtx4.columns)
mtx5_columns = list(mtx5.columns)
mtx6_columns = list(mtx6.columns)
mtx7_columns = list(mtx7.columns)
mtx8_columns = list(mtx8.columns)
mtx9_columns = list(mtx9.columns)

## Merge the variables on three datasets by "+" and remove duplicated values by setting "set"
columns_list = list(set(mtx_columns + mtx1_columns + mtx2_columns +
                        mtx3_columns + mtx4_columns + mtx5_columns +
                        mtx6_columns + mtx7_columns + mtx8_columns + mtx9_columns)) 

## Not to consider the "cell" value since it is working as key on the further manipulation
columns_list.remove("Cell") 
print(columns_list) 

## Run the loop for the getting average

result = pd.DataFrame(mtx["Cell"]) ## Build new data frame, and take the cell row.


for i in columns_list:
    count = 0 # Count the frequency of TFs in whole datasets
    is_in_mtx = False
    is_in_mtx1 = False
    is_in_mtx2 = False
    is_in_mtx3 = False
    is_in_mtx4 = False
    is_in_mtx5 = False
    is_in_mtx6 = False
    is_in_mtx7 = False
    is_in_mtx8 = False
    is_in_mtx9 = False

## create ith row and the default value is 0
    result[i] = 0 

## If variable i is included in the variable of mtx,
    if i in mtx.columns:
        count +=1
        result[i] = result[i] + mtx[i].astype("float") # ith row of mtx is changing as float format

## If variable i is included in the variable of mtx1,
    if i in mtx1.columns:
        count +=1
        result[i] = result[i] + mtx1[i].astype("float") # ith row of mtx is changing as float format

## If variable i is included in the variable of mtx2,
    if i in mtx2.columns:
        count +=1
        result[i] = result[i] + mtx2[i].astype("float")

    if i in mtx3.columns:
        count +=1
        result[i] = result[i] + mtx3[i].astype("float")

    if i in mtx4.columns:
        count +=1
        result[i] = result[i] + mtx4[i].astype("float")

    if i in mtx5.columns:
        count +=1
        result[i] = result[i] + mtx5[i].astype("float")

    if i in mtx6.columns:
        count +=1
        result[i] = result[i] + mtx6[i].astype("float")

    if i in mtx7.columns:
        count +=1
        result[i] = result[i] + mtx7[i].astype("float")

    if i in mtx8.columns:
        count +=1
        result[i] = result[i] + mtx8[i].astype("float")

    if i in mtx9.columns:
        count +=1
        result[i] = result[i] + mtx9[i].astype("float")        

    result[i] = result[i]/count # divide ith row of result by count to get the average

result.to_csv("AUC_Average.csv", index=False) # Export the result file as .csv file format.
