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



def GridSearch_RF(df, SAVE):
  import itertools

  #Make header array to save score from each random balanced replicate. Size = num_df
  results = np.array(['n_estimators', 'max_depth', 'max_features', 'criterion', 'accuracy', 'macro_f1'])
  
  # Set parameters to sweep in grid search
  n_estimators_list = [10, 50, 100, 500]
  max_depth_list = [2, 3, 5, 10, 20, 50, 100]
  max_features_list = [0.1, 0.25, 0.5, 0.75, .99, 'sqrt', 'log2']
  criterion_list = ['gini', 'entropy']
  
  for (n_estimators, max_depth, max_features, criterion) in itertools.product(n_estimators_list, max_depth_list, max_features_list, criterion_list):
    accuracy_l = np.array([])
    macro_f1_l = np.array([])

    for j in range(10):
        
      #Make balanced dataset with random negative examples drawn
      random_neg = all_neg.sample(n=pos_size, random_state = j)     #take random subset of negative examples - same size as positive
      df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe

      accuracy, macro_f1, importances = RandomForest(df, n_estimators, max_depth, max_features, criterion)

      # Add accuracy and f1 to results arrays
      accuracy_l = np.insert(accuracy_l, 0, accuracy)
      macro_f1_l = np.insert(macro_f1_l, 0, macro_f1) 

    new_row = [n_estimators, max_depth, max_features, criterion, np.mean(accuracy), float(np.mean(macro_f1))]
    results = np.vstack([results, new_row])

  r = pd.DataFrame(results[1:,:], columns = results[0,:])
  r = r.sort_values(['macro_f1','accuracy'], 0, ascending = [False, False])
  print("Top runs from the parameter sweep:")
  print(r.head(5))
  outName = SAVE + "_GridSearch"
  r.to_csv(outName, sep = "\t", index=False)

  return r['n_estimators'].iloc[0], r['max_depth'].iloc[0], r['max_features'].iloc[0], r['criterion'].iloc[0]
    
def GridSearch_LinearSVC(df, SAVE):
  import itertools

  #Make header array to save score from each random balanced replicate. Size = num_df
  results = np.array(['C', 'loss', 'max_iter', 'accuracy', 'macro_f1'])
  
  # Set parameters to sweep in grid search
  C_list = [0.01, 0.1, 0.5, 1, 10, 50, 100]
  loss_list = ['hinge', 'squared_hinge']
  max_iter_list=[10, 100,1000]
  
  for (C, loss, max_iter) in itertools.product(C_list, loss_list, max_iter_list):
    accuracy_l = np.array([])
    macro_f1_l = np.array([])

    for j in range(10):
        
      #Make balanced dataset with random negative examples drawn
      random_neg = all_neg.sample(n=pos_size, random_state = j)     #take random subset of negative examples - same size as positive
      df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe
      
      accuracy, macro_f1, importances = LinearSVC(df, C, loss, max_iter)
      
      # Add accuracy and f1 to results arrays
      accuracy_l = np.insert(accuracy_l, 0, accuracy)
      macro_f1_l = np.insert(macro_f1_l, 0, macro_f1) 

    new_row = [C, loss, max_iter, np.mean(accuracy), np.mean(macro_f1)]
    results = np.vstack([results, new_row])

  r = pd.DataFrame(results[1:,:], columns = results[0,:])
  r = r.sort_values(['accuracy', 'macro_f1'], 0, ascending = [False, False])
  print("Top runs from the parameter sweep:")
  print(r.head(5))
  outName = SAVE + "_GridSearch"
  r.to_csv(outName, sep = "\t", index=False)

  return r['C'].iloc[0], r['loss'].iloc[0], r['max_iter'].iloc[0]

def RandomForest(df, n_estimators, max_depth, max_features, criterion):
  from sklearn.ensemble import RandomForestClassifier
  from sklearn.model_selection import cross_val_predict
  from sklearn.metrics import accuracy_score, f1_score

  y = df['Class']   
  x = df.drop(['Class'], axis=1)  
  
  try:
    max_depth = int(max_depth)
  except:
    pass
  try:
    max_features = float(max_features)
  except:
    pass
  
  #Define the model
  clf = RandomForestClassifier(n_estimators=int(n_estimators), max_depth=max_depth, max_features=max_features, criterion=criterion)
  clf = clf.fit(x, y)
  #Obtain the predictions using 10 fold cross validation (uses KFold cv by default):
  cv_predictions = cross_val_predict(estimator=clf, X=x, y=y, cv=10)

  # Calculate the accuracy score and f1_score & add to 
  accuracy = accuracy_score(y, cv_predictions)
  macro_f1 = f1_score(y, cv_predictions)
  
  importances = clf.feature_importances_
  return accuracy, macro_f1, importances


def LinearSVC(df, C, loss, max_iter):
  from sklearn.svm import LinearSVC
  from sklearn.model_selection import cross_val_predict
  from sklearn.metrics import accuracy_score, f1_score

  y = df['Class']   
  x = df.drop(['Class'], axis=1)  
  
  #Define the model
  clf = LinearSVC(C=float(C), loss=loss, penalty='l2', max_iter=int(max_iter))
  clf = clf.fit(x, y)

  #Obtain the predictions using 10 fold cross validation (uses KFold cv by default):
  cv_predictions = cross_val_predict(estimator=clf, X=x, y=y, cv=10)

  # Calculate the accuracy score and f1_score & add to 
  accuracy = accuracy_score(y, cv_predictions)
  macro_f1 = f1_score(y, cv_predictions)
  
  importances = clf.coef_
  return accuracy, macro_f1, importances
 

if __name__ == "__main__":
    
  # Default code parameters
  neg, pos, n, FEAT, SAVE, GS, ALG = int(0), int(1), 50, 'all', 'test', 'False', 'RF' 

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



  ####### Run parameter sweep using a grid search #######
  if GS == 'True' or GS == 'T' or GS == 'y' or GS == 'yes':
    imp = "False"
    if ALG == "RF":
      n_estimators, max_depth, max_features, criterion = GridSearch_RF(df, SAVE)
    elif ALG == "SVC":
      C, loss, max_iter = GridSearch_LinearSVC(df, SAVE)


  ####### ML Pipeline #######

  acc = np.array([])
  f1 = np.array([])
  imp_array = list(df.columns.values)[1:]

  for j in range(n):
    
    #Make balanced dataset with random negative examples drawn
    random_neg = all_neg.sample(n=pos_size, random_state = j)     #take random subset of negative examples - same size as positive
    df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe
      
    if ALG == "RF":
      accuracy, macro_f1, importances = RandomForest(df, n_estimators, max_depth, max_features, criterion)
    elif ALG == "SVC":
      accuracy, macro_f1, importances = LinearSVC(df, C, loss, max_iter)
    
    # Add accuracy, f1, and importance scores to results arrays
    acc = np.insert(acc, 0, accuracy)
    f1 = np.insert(f1, 0, macro_f1)
    imp_array = np.vstack((imp_array, importances))

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
