#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import dependencies
import numpy as np
import pandas as pd
#from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score as acc
from mlxtend.feature_selection import SequentialFeatureSelector as sfs


# In[ ]:


# Read data
df = pd.read_csv('~/Documents/PhD_Project/Immune_Cells/Fragments_Annotated/Mono_NG', sep='\t')
print(df)


# In[ ]:


# Format data
df = df[['Propagation', 'eQTL_frequency', 'CAGEseq_frequency', 'Active']]
print(df)
df['Bin'] = [1 if x >= 0.5 else 0 for x in df['Active']]
print(df.loc[df['Active'] >= 0.5])
print(df.loc[df['Active'] <= 0.5])
df = df[['Propagation', 'eQTL_frequency', 'CAGEseq_frequency', 'Bin']]
print(df)


# In[ ]:


# Train/test split
# TEST DATA df = pd.read_csv('/home/maninder/Downloads/winequality-white.csv', sep = ';')

X_train, X_test, y_train, y_test = train_test_split(
    df.values[:,:-1],
    df.values[:,-1:],
    test_size=0.25,
    random_state=42)

y_train = y_train.ravel()
y_test = y_test.ravel()

print('Training dataset shape:', X_train.shape, y_train.shape)
print('Testing dataset shape:', X_test.shape, y_test.shape)
#print(X_train)
#print(X_test)
print(y_train)
print(y_test)
#print(df.values)


# In[ ]:


# Build RF classifier to use in feature selection
clf = RandomForestClassifier(n_estimators=100, n_jobs=-1)

# Build step forward feature selection
sfs1 = sfs(clf,
           k_features=3,
           forward=True,
           floating=False,
           verbose=2,
           scoring='accuracy',
           cv=5)

# Perform SFFS
sfs1 = sfs1.fit(X_train, y_train)


# In[ ]:


# Which features?
feat_cols = list(sfs1.k_feature_idx_)
print(feat_cols)


# In[ ]:


clf = RandomForestClassifier(n_estimators=1000, random_state=42, max_depth=4)
clf.fit(X_train[:, feat_cols], y_train)

y_train_pred = clf.predict(X_train[:, feat_cols])
print('Training accuracy on selected features: %.3f' % acc(y_train, y_train_pred))

y_test_pred = clf.predict(X_test[:, feat_cols])
print('Testing accuracy on selected features: %.3f' % acc(y_test, y_test_pred))


# In[104]:


# Compare with an SVM
from sklearn.svm import SVC
svm = SVC(kernel="linear", C=0.025, random_state=101)
svm.fit(X_train, y_train)
y_pred=svm.predict(X_test)


# In[107]:


print(' acc '% acc(y_pred))

