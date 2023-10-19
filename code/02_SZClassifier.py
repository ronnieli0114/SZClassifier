# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 16:18:29 2023
@author: ronnieli
"""

import pandas as pd
import numpy as np
import os
import pyarrow.feather as feather
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score

class SZClassifier():
    
    def __init__(self, data_dir, n_top_cpgs, direction='all'):
        self.data_dir = data_dir
        self.n_top_cpgs = n_top_cpgs
        self.direction = direction
     
    def select_cpgs(self):
        
        ttest = pd.read_csv(os.path.join(self.data_dir,'sz_classifier','SZ_vs_control_ttest.txt'), sep='\t', header=0)
        if self.direction == 'all':
            print('Selecting all probes')
        elif self.direction == 'hyper':
            print('Selecting HYPERmethylated probes only')
            ttest = ttest[ttest['mean_diff'] > 0]
        elif self.direction == 'hypo':
            print('Selecting HYPOmethylated probes only')
            ttest = ttest[ttest['mean_diff'] < 0]
        
        ttest.sort_values('p_adj', ascending=True, inplace=True)
        
        if self.n_top_cpgs == 'all':
            cpg_list = ttest['cpg'].tolist()
        else:
            cpg_list = ttest.iloc[:self.n_top_cpgs,:]['cpg'].tolist()

        return cpg_list
        
    def load_Xy(self): 
        
        cpg_list = self.select_cpgs()
        
        print(f'Selecting {len(cpg_list)} top CpGs from a t-test...')
        cpg_meth = feather.read_feather(os.path.join(self.data_dir,'sz_classifier','GSE147221_methylation.ftr'))
        subj = pd.read_csv(os.path.join(self.data_dir,'sz_classifier','GSE147221_subject_info.txt'), sep='\t', header=0)
        print('Loaded all raw data')
        
        cpg_meth.set_index('cpg', inplace=True)
        features = cpg_meth[cpg_meth.index.isin(cpg_list)].T    
        labels = subj.loc[:,['sample_id','target']].set_index('sample_id')
        
        Xy = features.merge(labels, left_index=True, right_index=True)
        Xy['target'] = Xy['target'].apply(lambda s: s.upper().replace('CASE','1').replace('CONTROL','0'))
        Xy['target'] = Xy['target'].astype(int)
    
        X = Xy.drop('target', axis=1)
        y = Xy['target']
        
        return X, y
        
    def classify(self):
        
        result_dict = {}
        
        X, y = self.load_Xy()
        
        X = X.values
        y = y.values.reshape(-1,)
        
        skf = StratifiedKFold(n_splits=5)
        clf = RandomForestClassifier(n_estimators=150, max_depth=4, random_state=2023)
        # clf = SVC(C=1.0, kernel='rbf', gamma='auto', probability=True)
        scaler = StandardScaler()
        
        f1_list = []
        auc_list = []
        mcc_list = []
        
        for train_idx, test_idx in skf.split(X, y):
            
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
            
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            
            f1 = f1_score(y_test, y_pred, average='macro')
            auc = roc_auc_score(y_test, y_pred)
            mcc = matthews_corrcoef(y_test, y_pred)
            
            f1_list.append(f1)
            auc_list.append(auc)
            mcc_list.append(mcc)
        
        f1_avg = np.around(np.mean(f1_list), 3)
        auc_avg = np.around(np.mean(auc_list), 3)
        mcc_avg = np.around(np.mean(mcc_list), 3)
        
        result_dict['F1'] = f1_avg
        result_dict['AUC'] = auc_avg
        result_dict['MCC'] = mcc_avg
        
        return pd.Series(result_dict)
    
    def classify_with_feature_selection(self, K):
        
        X, y = self.load_Xy()
        
        # select K best features
        print(f'Selecting {K} best features...')
        kbest = SelectKBest(f_classif, k=K)
        X_trans = kbest.fit_transform(X, y)
        selected_features = kbest.get_feature_names_out()

        # classify
        result_dict = {}
        
        skf = StratifiedKFold(n_splits=5)
        clf = SVC(C=1.0, kernel='rbf', gamma='auto', probability=True)
        scaler = StandardScaler()
        
        f1_list = []
        auc_list = []
        mcc_list = []
        
        for train_idx, test_idx in skf.split(X_trans, y):
            
            X_train, X_test = X_trans[train_idx], X_trans[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
            
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            
            f1 = f1_score(y_test, y_pred, average='macro')
            auc = roc_auc_score(y_test, y_pred)
            mcc = matthews_corrcoef(y_test, y_pred)
            
            f1_list.append(f1)
            auc_list.append(auc)
            mcc_list.append(mcc)
        
        f1_avg = np.around(np.mean(f1_list), 3)
        auc_avg = np.around(np.mean(auc_list), 3)
        mcc_avg = np.around(np.mean(mcc_list), 3)
        
        result_dict['F1'] = f1_avg
        result_dict['AUC'] = auc_avg
        result_dict['MCC'] = mcc_avg
        
        return list(selected_features), pd.Series(result_dict)

model = SZClassifier(data_dir=r'E:/lab_data', n_top_cpgs=50, direction='all')   
feats, result = model.classify_with_feature_selection(K=10)
print(result.head())
