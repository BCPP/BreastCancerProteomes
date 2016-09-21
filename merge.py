import matplotlib.pyplot as plt
from matplotlib.mlab import PCA as mlabPCA
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

#print processed_numerical[0]


## start PCA algorithm
mlab_pca = mlabPCA(processed_numerical) 

fracs = mlab_pca.fracs
vector = mlab_pca.Y

## draw Proportion Fractions Plot
a=80
xc = []
for i in range(a):
 xc.append(i)


plt.plot(xc, fracs, 'o', markersize=7, color='blue', alpha=0.5, label='fracs')

plt.xlim([0,80])
plt.ylim([0,0.15])
plt.xlabel('Columns')
plt.ylabel('Proportion fraction')
plt.legend()
plt.title('The proportion of variance of each column')

plt.show()

#print('Fracs:\n')
#print fracs

## tuning point, choose 0.02, remove low porportion columns
select = 0
for f in fracs:
  if f > 0.02:  
     select += 1

afterPCA_numerical = mlab_pca.Y[0 : 12552, 0 : select]
#print processed_numerical[0]
#print afterPCA_numerical[0]
