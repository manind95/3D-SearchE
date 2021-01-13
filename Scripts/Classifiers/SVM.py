#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn import metrics


# In[ ]:


# Read data
df = pd.read_csv('~/Documents/PhD_Project/Immune_Cells/Fragments_Annotated/Mono_NG', sep='\t')


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


# Split the data into train and test where x is the features and y is the labels
x_train, x_test, y_train, y_test = train_test_split(
    df.values[:,:-1],
    df.values[:,-1:],
    test_size=0.25,
    random_state=42)

y_train = y_train.ravel()
y_test = y_test.ravel()

print('Training dataset shape:', x_train.shape, y_train.shape)
print('Testing dataset shape:', x_test.shape, y_test.shape)


# In[ ]:


# Generate the model (Specify the kernal = ('linear', 'poly', 'rbf'))
clf = svm.SVC(kernel = 'linear')
clf.fit(x_train, y_train)


# In[ ]:


# Generate the predictions
y_pred = clf.predict(x_test)


# In[ ]:


# Model Accuracy: how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))


# In[ ]:


# Model Precision: what percentage of positive tuples are labeled as such?
print("Precision:",metrics.precision_score(y_test, y_pred))


# In[ ]:


# Model Recall: what percentage of positive tuples are labelled as such?
print("Recall:",metrics.recall_score(y_test, y_pred))


# In[ ]:


# Model ROC curve
fpr, tpr, thresholds = "ROC:", metrics.roc_curve(y_test, y_pred, pos_label=None, sample_weight=None, drop_intermediate=True)
#roc_auc = metrics.auc(fpr, tpr)

# Plot: ggplot
#from ggplot import *
#df = pd.DataFrame(dict(fpr = fpr, tpr = tpr))
#ggplot(df, aes(x = 'fpr', y = 'tpr')) + geom_line() + geom_abline(linetype = 'dashed')


# In[ ]:




