"""
PURPOSE:
Bianary classifications using RandomForestClassifier and LinearSVC implemented in sci-kit learn. 

Optional grid seach function (-gs True) to run a parameter sweep on a subset of the balanced datasets

To access pandas, numpy, and sklearn packages on HPS first run:
$ export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH


INPUTS:
  
  REQUIRED:
  -df       Feature & class dataframe for ML. See "" for an example dataframe
  -save     Unique save name (caution - will overwrite!)
  
  OPTIONAL:
  -feat     Import file with list of features to keep if desired. Default: keep all.
  -gs       Set to True if parameter sweep is desired. Default = False
  -alg      Algorithm to use. Currently available: RandomForest (RF)(Default), LinearSVC (SVC)
  -n        Number of random balanced datasets to run. Default = 50
  -pos      String for what codes for the positive example (i.e. UUN) Default = 1
  -neg      String for what codes for the negative example (i.e. NNN) Default = 0


OUTPUT:
  -SAVE_imp           Importance scores for each feature
  -SAVE_GridSearch    Results from parameter sweep sorted by F1
  -RESULTS_XX.txt     Accumulates results from all ML runs done in a specific folder - use unique save names! XX = RF or SVC

"""
import sys
import pandas as pd
import numpy as np
from math import sqrt
from sklearn.metrics import (precision_recall_curve, f1_score)
import time


def GridSearch_RF(df):
  from sklearn.model_selection import GridSearchCV
  from sklearn.ensemble import RandomForestClassifier
  start_time = time.time()

  parameters = {'n_estimators':(50, 100, 500), 'max_depth':(2, 3, 5, 10, 50, 100), 'max_features': (0.1, 0.25, 0.5, 0.75, 0.9999, 'sqrt', 'log2')}
  #parameters = {'n_estimators':(10, 50), 'max_depth':(2, 5), 'max_features': (0.1, 0.5)}
  
  for j in range(10):
    # Build random balanced dataset and define x & y
    random_neg = all_neg.sample(n=pos_size, random_state = j)     #take random subset of negative examples - same size as positive
    df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe
    
    y = df['Class']   
    x = df.drop(['Class'], axis=1) 

    # Build model, run grid search with 10-fold cross validation and fit
    rf = RandomForestClassifier()
    clf = GridSearchCV(rf, parameters, cv = 10, n_jobs = n_jobs, pre_dispatch=2*n_jobs)
    clf.fit(x, y)
    j_results = pd.DataFrame(clf.cv_results_)
    
    if j == 0:
      results = j_results[['param_max_depth', 'param_max_features', 'param_n_estimators','mean_test_score']]
    else:
      results_temp = j_results[['param_max_depth', 'param_max_features', 'param_n_estimators','mean_test_score']]
      results = pd.merge(results, results_temp, on=['param_max_depth', 'param_max_features', 'param_n_estimators'])
  
  # Calculate average test score and sort
  col_list= [col for col in results.columns if 'mean_test_score' in col]
  results['average_mean_test_score'] = results[col_list].mean(axis=1)
  results = results.sort_values(['average_mean_test_score'], 0, ascending = False)

  print("Parameter sweep time: %f seconds" % (time.time() - start_time))
  outName = SAVE + "_GridSearch"
  results.to_csv(outName, sep = "\t", index=False, columns = ['param_max_depth', 'param_max_features', 'param_n_estimators', 'average_mean_test_score'])

  return results['param_n_estimators'].iloc[0], results['param_max_depth'].iloc[0], results['param_max_features'].iloc[0]


def GridSearch_LinearSVC(df):
  from sklearn.model_selection import GridSearchCV
  from sklearn.svm import LinearSVC
  start_time = time.time()

  parameters = {'C':(0.01, 0.1, 0.5, 1, 10, 50, 100), 'loss':('hinge', 'squared_hinge'), 'max_iter':(10,100,1000)}
  
  for j in range(20):
    # Build random balanced dataset and define x & y
    random_neg = all_neg.sample(n=pos_size, random_state = j)     #take random subset of negative examples - same size as positive
    df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe
    
    y = df['Class']   
    x = df.drop(['Class'], axis=1) 

    # Build model, run grid search with 10-fold cross validation and fit
    svc = LinearSVC()
    clf = GridSearchCV(svc, parameters, cv = 10, n_jobs = 100)
    clf.fit(x, y)
    j_results = pd.DataFrame(clf.cv_results_)
    
    if j == 0:
      results = j_results[['param_C', 'param_loss', 'param_max_iter','mean_test_score']]
    else:
      results_temp = j_results[['param_C', 'param_loss', 'param_max_iter','mean_test_score']]
      results = pd.merge(results, results_temp, on=['param_C', 'param_loss', 'param_max_iter'])
  
  # Calculate average test score and sort
  col_list= [col for col in results.columns if 'mean_test_score' in col]
  results['average_mean_test_score'] = results[col_list].mean(axis=1)
  results = results.sort_values(['average_mean_test_score'], 0, ascending = False)

  print("Parameter sweep time: %f seconds" % (time.time() - start_time))
  outName = SAVE + "_GridSearch"
  results.to_csv(outName, sep = "\t", index=False, columns = ['param_C', 'param_loss', 'param_max_iter', 'average_mean_test_score'])

  return results['param_C'].iloc[0], results['param_loss'].iloc[0], results['param_max_iter'].iloc[0]


def RandomForest(df, n_estimators, max_depth, max_features, criterion, n_jobs):
  from sklearn.ensemble import RandomForestClassifier
  from sklearn.model_selection import cross_val_predict
  from sklearn.metrics import accuracy_score, f1_score

  y = df['Class']   
  x = df.drop(['Class'], axis=1) 
  
  #Define the model
  clf = RandomForestClassifier(n_estimators=int(n_estimators), max_depth=max_depth, max_features=max_features, criterion=criterion, n_jobs=n_jobs)
  clf = clf.fit(x, y)
  
  #Obtain the predictions using 10 fold cross validation (uses KFold cv by default):
  cv_predictions = cross_val_predict(estimator=clf, X=x, y=y, cv=10)
  #probability_pos_clf = clf.predict_proba()

  # Calculate the accuracy score and f1_score & add to 
  accuracy = accuracy_score(y, cv_predictions)
  macro_f1 = f1_score(y, cv_predictions)
  
  importances = clf.feature_importances_
  return accuracy, macro_f1, importances #, cv_predictions


def LinearSVC(df, C, loss, max_iter, n_jobs):
  from sklearn.svm import LinearSVC
  from sklearn.model_selection import cross_val_predict
  from sklearn.metrics import accuracy_score, f1_score

  y = df['Class']   
  x = df.drop(['Class'], axis=1)  
  
  #Define the model
  clf = LinearSVC(C=float(C), loss=loss, penalty='l2', max_iter=int(max_iter))
  clf = clf.fit(x, y)

  #Obtain the predictions using 10 fold cross validation (uses KFold cv by default):
  cv_predictions = cross_val_predict(estimator=clf, X=x, y=y, cv=10, n_jobs = n_jobs)

  # Calculate the accuracy score and f1_score & add to 
  accuracy = accuracy_score(y, cv_predictions)
  macro_f1 = f1_score(y, cv_predictions)
  
  importances = clf.coef_
  return accuracy, macro_f1, importances
 

def PR_Curve(y_pred, SAVE):
  import matplotlib.pyplot as plt
  plt.switch_backend('agg')

  
  f1_summary = f1_score(y_true, y_pred_round)
  plt.plot(recall, precision, color='navy', label='Precision-Recall curve')
  plt.xlabel('Recall')
  plt.ylabel('Precision')
  plt.ylim([0.0, 1.05])
  plt.xlim([0.0, 1.0])
  plt.title('PR Curve: %s\nSummary F1 = %0.2f' % (SAVE, f1_summary))
  plt.show()
  # save a PDF file named for the CSV file (but in the current directory)
  filename = SAVE + "_PRcurve.png"
  plt.savefig(filename)



if __name__ == "__main__":
    
  # Default code parameters
  neg, pos, n, FEAT, SAVE, GS, ALG, PR, n_jobs = int(0), int(1), 50, 'all', 'test', 'False', 'RF', "False", 100

  # Default Random Forest parameters
  n_estimators, max_depth, max_features, criterion = 500, 10, "sqrt", "gini"

  # Default Linear SVC parameters
  C, loss, max_iter = 1, 'hinge', "500"

  for i in range (1,len(sys.argv),2):
        if sys.argv[i] == "-df":
          DF = sys.argv[i+1]
        if sys.argv[i] == '-save':
          SAVE = sys.argv[i+1]
        if sys.argv[i] == '-feat':
          FEAT = sys.argv[i+1]
        if sys.argv[i] == "-gs":
          GS = sys.argv[i+1]
        if sys.argv[i] == '-neg':
          neg = sys.argv[i+1]
        if sys.argv[i] == "-pos":
          pos = sys.argv[i+1]
        if sys.argv[i] == "-n":
          n = int(sys.argv[i+1])
        if sys.argv[i] == "-alg":
          ALG = sys.argv[i+1]
        if sys.argv[i] == "-criterion":
          criterion = sys.argv[i+1]
        if sys.argv[i] == "-PR":
          PR = sys.argv[i+1]
        if sys.argv[i] == "-n_jobs":
          n_jobs = int(sys.argv[i+1])


  if len(sys.argv) <= 1:
    print(__doc__)
    exit()
  

  ####### Load Dataframe & Pre-process #######

  if isinstance(DF, str):
    df = pd.read_csv(DF, sep='\t', index_col = 0)
  else:
    df = DF

  
  # If list of features to include in analysis given (-feat), filter out all other features
  if FEAT != 'all':
    with open(FEAT) as f:
      features = f.read().splitlines()
      features = ['Class'] + features
    df = df.loc[:,features]

  # Recode class as 1 for positive and 0 for negative, then divide into two dataframes.
  df["Class"] = df["Class"].replace(pos, 1)
  df["Class"] = df["Class"].replace(neg, 0)

  # Determine number of positive and negative instances
  n_features = len(list(df.columns.values)[1:])
  all_pos = df[df.Class == 1]
  pos_size = all_pos.shape[0]
  all_neg = df[df.Class == 0]
  print("%d positive examples and %d negative examples in dataframe" % (pos_size, all_neg.shape[0]))

  SAVE = SAVE + "_" + ALG

  ####### Run parameter sweep using a grid search #######
  if GS == 'True' or GS == 'T' or GS == 'y' or GS == 'yes':
    imp = "False"
    if ALG == "RF":
      n_estimators, max_depth, max_features = GridSearch_RF(df)
      print("Parameters selected: n_estimators=%s, max_depth=%s, max_features=%s" % (str(n_estimators), str(max_depth), str(max_features)))
    elif ALG == "SVC":
      C, loss, max_iter = GridSearch_LinearSVC(df)
      print("Parameters selected: C=%s, loss=%s, max_iter=%s" % (str(C), str(loss), str(max_iter)))


  ####### ML Pipeline #######
  start_time = time.time()
  acc = np.array([])
  f1 = np.array([])
  positives, negatives = np.ones(pos_size), np.zeros(pos_size)
  y_pred = y_true = np.hstack((positives, negatives))
  imp_array = list(df.columns.values)[1:]

  for j in range(n):
    
    #Make balanced dataset with random negative examples drawn
    random_neg = all_neg.sample(n=pos_size, random_state = j)     #take random subset of negative examples - same size as positive
    df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe

    if ALG == "RF":
      accuracy, macro_f1, importances = RandomForest(df, n_estimators, max_depth, max_features, criterion, n_jobs)
    elif ALG == "SVC":
      accuracy, macro_f1, importances = LinearSVC(df, C, loss, max_iter, n_jobs)
  
    # Calculate precision and recall for the run
    #precision, recall, thresholds = precision_recall_curve(y_true, cv_predictions, pos_label=1)
    #print(precision, recall, thresholds)

    # Add accuracy, f1, and importance scores to results arrays
    acc = np.insert(acc, 0, accuracy)
    f1 = np.insert(f1, 0, macro_f1)
    imp_array = np.vstack((imp_array, importances))
      
    
    #y_pred = np.vstack((y_pred, cv_predictions))
  # Make PR Curve (also code to ID which instances are being predicted)
  print("ML Pipeline time: %f seconds" % (time.time() - start_time))

  if PR == "True" or PR == "T":
    PR_Curve(y_pred, SAVE)

  # Calculate mean importance scores and sort
  imp_df = pd.DataFrame(imp_array[1:,:], columns=imp_array[0,:], dtype = float)  
  imp_mean_df = pd.DataFrame(imp_df.mean(axis=0, numeric_only=True), columns = ['importance'], index = imp_df.columns.values)
  imp_mean_df = imp_mean_df.sort_values('importance', 0, ascending = False)
  print("Top five most important features:")
  print(imp_mean_df.head(5))
  
  imp_out = SAVE + "_imp"
  imp_mean_df.to_csv(imp_out, sep = "\t", index=True)

  # Calculate accuracy and f1 stats
  acc_stdv = np.std(acc)
  acc_SE = np.std(acc)/sqrt(n)
  f1_stdv = np.std(f1)
  f1_SE = np.std(f1)/sqrt(n)
  
  print("\nML Results: \nAccuracy: %03f (+/- stdev %03f)\nF1 (macro): %03f (+/- stdev %03f)" % (np.mean(acc), acc_stdv, np.mean(f1), f1_stdv))
   
  if ALG == "RF":
    open("RESULTS_RF.txt",'a').write('%s\t%i\t%i\t%i\t%i\t%s\t%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n' % (SAVE, pos_size, all_neg.shape[0], n_features, n, str(n_estimators), str(max_depth), str(max_features), str(criterion), np.mean(acc), acc_stdv, acc_SE, np.mean(f1), f1_stdv, f1_SE))
    print ('\nColumn Names in RF_RESULTS.txt output: Run_Name, #Pos_Examples, #Neg_Examples, #Features, #Random_df_Reps, n_estimators, max_depth, max_features, criterion, Accuracy, Accuracy_Stdev, Accuracy_SE, F1, F1_Stdev, F1_SE')

  elif ALG == "SVC":
    open("RESULTS_SVC.txt",'a').write('%s\t%i\t%i\t%i\t%i\t%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n' % (SAVE, pos_size, all_neg.shape[0], n_features, n, str(C), str(loss), str(max_iter), np.mean(acc), acc_stdv, acc_SE, np.mean(f1), f1_stdv, f1_SE))
    print ('\nColumn Names in LinearSVC_RESULTS.txt output: Run_Name, #Pos_Examples, #Neg_Examples, #Features, #Random_df_Reps, C, loss, max_iter, Accuracy, Accuracy_Stdev, Accuracy_SE, F1, F1_Stdev, F1_SE')
