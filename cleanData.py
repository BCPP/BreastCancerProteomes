
# coding: utf-8

# In[91]:

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
from sklearn import metrics
import re
from sklearn.preprocessing import Imputer
import numpy as np
from numpy import random


#set dataset path
dataset_path = "./77_cancer_proteomes_CPTAC_itraq.csv"
clinical_info = "./clinical_data_breast_cancer.csv"
pam50_proteins = "./PAM50_proteins.csv"

## Load data
data = pd.read_csv(dataset_path,header=0,index_col=0)
clinical = pd.read_csv(clinical_info,header=0,index_col=0)## holds clinical information about each patient/sample
pam50 = pd.read_csv(pam50_proteins,header=0)

## Drop unused information columns
data.drop(['gene_symbol','gene_name'],axis=1,inplace=True)

## Drop three healthy samples
data.drop(['263d3f-I.CPTAC','blcdb9-I.CPTAC','c4155b-C.CPTAC'], axis = 1, inplace = True);
 
## Change the protein data sample names to a format matching the clinical data set
data.rename(columns=lambda x: "TCGA-%s" % (re.split('[_|-|.]',x)[0]) if bool(re.search("TCGA",x)) is True else x,inplace=True)

## Numerical data for the algorithm
data.loc[:,[x for x in data.columns if bool(re.search("TCGA",x)) == True]]

## Impute missing values by median of the columns 
imputer = Imputer(missing_values='NaN', strategy='median', axis=0)
imputer = imputer.fit(data)
processed_numerical = imputer.transform(data)

#print processed_numerical.shape

#save cleaned data as numpy array to csv
#np.savetxt(
#    'cleaned_data.csv',           # file name
#    processed_numerical,      # array to save
#    fmt='%f',             # formatting, 2 digits in this case
#    delimiter=',',          # column delimiter
#    newline='\n',           # new line character
#    comments='# ',          # character to use for comments
#    header='Cleaned Data generated by numpy') 





