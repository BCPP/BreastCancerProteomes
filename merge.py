import matplotlib.pyplot as plt
from matplotlib.mlab import PCA as mlabPCA
import pandas as pd
from sklearn.cluster import KMeans
from sklearn import metrics
import re
from sklearn.preprocessing import Imputer
import numpy as np
from numpy import random
import seaborn as sns

#set dataset path
dataset_path = "./77_cancer_proteomes_CPTAC_itraq.csv"
clinical_info = "./clinical_data_breast_cancer.csv"
pam50_proteins = "./PAM50_proteins.csv"

## Load data
data = pd.read_csv(dataset_path,header=0,index_col=0)
clinical = pd.read_csv(clinical_info,header=0,index_col=0)

## holds clinical information about each patient/sample
pam50 = pd.read_csv(pam50_proteins,header=0)

## Drop unused information columns
data.drop(['gene_symbol','gene_name'],axis=1,inplace=True)

## Drop three healthy samples
data.drop(['263d3f-I.CPTAC','blcdb9-I.CPTAC','c4155b-C.CPTAC'], axis = 1, inplace = True);
 
## Change the protein data sample names to a format matching the clinical data set
data.rename(columns=lambda x: "TCGA-%s" % (re.split('[_|-|.]',x)[0]) if bool(re.search("TCGA",x)) is True else x,inplace=True)

data = data.transpose()
clinical = clinical.loc[[x for x in clinical.index.tolist() if x in data.index],:]

#merge 77 samples with clinical data
data_for_merge = data
merged = data_for_merge.merge(clinical,left_index=True,right_index=True)

## Numerical data for the algorithm
#data.loc[:,[x for x in data.columns if bool(re.search("TCGA",x)) == True]]

merged.loc[:,[x for x in data.columns if bool(re.search("NP_|XP_",x)) == True]]
processed_merged = merged.ix[:,merged.columns.isin(pam50['RefSeqProteinID'])]

## Impute missing values by median of the rows
imputer = Imputer(missing_values='NaN', strategy='median', axis=1)
imputer = imputer.fit(processed_merged)
processed_merged = imputer.transform(processed_merged)


## start PCA algorithm
mlab_pca = mlabPCA(processed_merged) 

fracs = mlab_pca.fracs
#print fracs.shape



vector = mlab_pca.Y

## draw Proportion Fractions Plot
a = 43
xc = []
for i in range(a):
 xc.append(i)

plt.plot(xc, fracs, 'o', markersize=7, color='blue', alpha=0.5, label='fracs')

plt.xlim([0,43])
plt.ylim([0,0.3])
plt.xlabel('Proteomes')
plt.ylabel('Proportion fraction')
plt.legend()
plt.title('The proportion of variance of each column')

#plt.show()



#print fracs

## tuning point, choose 0.02, remove low porportion columns
select = 0
for f in fracs:
  if f > 0.01:  
     select += 1

afterPCA_numerical = mlab_pca.Y[0 : 80, 0 : select]

#print afterPCA_numerical.shape

n_clusters = [2,3,4,5,6,7,8,10]

def compare_k_means(k_list,data):
    ## Run clustering with different k and check the metrics
    for k in k_list:
        clusterer = KMeans(n_clusters=k, n_jobs=4)
        clusterer.fit(data)
        ## The higher (up to 1) the better
        print("Silhouette Coefficient for k == %s: %s" % (
        k, round(metrics.silhouette_score(data, clusterer.labels_), 4)))
        ## The higher (up to 1) the better
        print("Homogeneity score for k == %s: %s" % (
        k, round(metrics.homogeneity_score(merged['PAM50 mRNA'], clusterer.labels_),4)))
        print("------------------------")

#compare_k_means(n_clusters, processed_numerical)
#compare_k_means(n_clusters, processed_merged)
compare_k_means(n_clusters, afterPCA_numerical)


clusterer_final = KMeans(n_clusters = 4, n_jobs = 4)
clusterer_final = clusterer_final.fit(afterPCA_numerical)
afterPCA_plot = pd.DataFrame(afterPCA_numerical)
afterPCA_plot['KMeans_cluster'] = clusterer_final.labels_
afterPCA_plot.sort_values('KMeans_cluster',axis = 0,inplace=False)
afterPCA_plot.index.name = 'Patient'
sns.heatmap(afterPCA_plot)
plt.title('KMeans_heatmap')
#plt.show()
