#!/usr/bin/python -u
from __future__ import division
import sys, getopt
import re
import numpy as np
from scipy import stats
import argparse, os, sys
import patsy
import math
import pandas as pd
import statsmodels.api as sm
import sklearn
from sklearn.model_selection import KFold # import KFold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeavePOut
from functools import reduce

k=sys.argv[1]
kf = KFold(n_splits=k)

###########################################################################
r=pd.read_csv("data/final_admixture_results_merged_egfr_kras.csv", sep="\t")
markers=r.columns

P={}
C={}
Z={}
O1={}
O2={}

for l in range(1, k+1):
    O1[l]=open("admixture_stats_k"+str(l)+"_all.txt", 'w')
    O2[l]=open("admixture_results_k"+str(l)+"_all.txt", 'w')

l=0

##############################################################################
############### start to use from here ################
##############################################################################

for train_index, test_index in kf.split(r["EGFR"].values):
#train, test, y_train, y_test = train_test_split(r, r["EGFR"].values, test_size=0.4)
        train=r.iloc[train_index]
        test=r.iloc[test_index]

        train=train.dropna()
        test=test.dropna()
        test["risk_sum0"]=0
        test["risk_sum1"]=0
        test["risk_sum2"]=0
        test["risk_sum3"]=0

        print (len(train), len(test))
        print (len(train[train["EGFR"]==1]), len(test[test["EGFR"]==1]))
        for i in markers[1:-1]:
                if ":" in i:
                        print (i)

                        ############################################
                        ## discovery ancestry associated loci ######
                        ############################################
                        new=train[[i, "EGFR", "AMR_mean", "cohort"]]
                        new.columns=["marker", "EGFR", "AMR_mean", "cohort"]
                        new["marker"] = new["marker"].astype(int)
                        #null_columns=new.columns[new.isnull().any()]
                        #print new[new.isnull().any(axis=1)][null_columns].head()
                        f = "EGFR ~ marker + AMR_mean + cohort"
                        y, X = patsy.dmatrices(f, new, return_type='dataframe')
                        result = sm.Logit(y, X).fit()
                        P[i] = result.pvalues[2:3].values[0]
                        C[i] = result.params[2:3].values[0]
                        Z[i] = result.summary2().tables[1]['z'][2:3].values[0]

                        #########################################################
                        ## calculate polygeneic risk score for the test cohort ##
                        #########################################################
                        if (float(P[i])<0.00005) and (float(C[i])>0):
                                test["risk_sum0"]=test["risk_sum0"]+test[i].astype(int)*Z[i].astype(int)
                        if (float(P[i])<0.0001) and (float(C[i])>0):
                                test["risk_sum1"]=test["risk_sum1"]+test[i].astype(int)*Z[i].astype(int)

                        if (float(P[i])<0.001) and (float(C[i])>0):
                                test["risk_sum2"]=test["risk_sum2"]+test[i].astype(int)*Z[i].astype(int)

                        if (float(P[i])<0.05) and (float(C[i])>0):
                                test["risk_sum3"]=test["risk_sum3"]+test[i].astype(int)*Z[i].astype(int)

        l=l+1
        test[["#sample", "EGFR", "risk_sum0", "risk_sum1", "risk_sum2", "risk_sum3", "AMR_mean", "cohort"]].to_csv(O2[l], sep="\t")

        f = "EGFR ~ risk_sum1 + AMR_mean + cohort"
        y, X = patsy.dmatrices(f, test, return_type='dataframe')
        result = sm.Logit(y, X).fit()
        print (result.pvalues, file=O1[l])
        print (result.summary(), file=O1[l])

        f = "EGFR ~ risk_sum2 + AMR_mean + cohort"
        y, X = patsy.dmatrices(f, test, return_type='dataframe')
        result = sm.Logit(y, X).fit()
        print (result.pvalues, file=O1[l])
        print (result.summary(), file=O1[l])

        f = "EGFR ~ risk_sum3 + AMR_mean + cohort"
        y, X = patsy.dmatrices(f, test, return_type='dataframe')
        result = sm.Logit(y, X).fit()
        print (result.pvalues, file=O1[l])
        print (result.summary(), file=O1[l])

        f = "EGFR ~ risk_sum0 + AMR_mean + cohort"
        y, X = patsy.dmatrices(f, test, return_type='dataframe')
        result = sm.Logit(y, X).fit()
        print (result.pvalues, file=O1[l])
        print (result.summary(), file=O1[l])
