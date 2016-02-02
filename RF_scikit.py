"""
PURPOSE:
Run classification RF from sci-kit learn on a given dataframe
*When submitting jobs ask for 8 nodes! 

INPUT:
  -df       FASTA file with positive examples
  -save     Save name (will overwrite some results files if the same as other names in the directory youre exporting to)

OUTPUT:
  -SAVE_RF_results.txt    Results from RF runs
  -RESULTS.txt            Final results get added to this file: Run Name, # Features, # Reps (different Neg Datasets), CV, F_measure, StDev, SE
"""

import pandas as pd
import numpy as np
import sys
from math import sqrt


for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-df":
        DF = sys.argv[i+1]
      if sys.argv[i] == '-save':
        SAVE = sys.argv[i+1]
      if sys.argv[i] == '-neg':
        neg = sys.argv[i+1]
      if sys.argv[i] == "-pos":
        pos = sys.argv[i+1]


def RandomForest(DF, SAVE, pos, neg):

  from sklearn import cross_validation
  from sklearn.preprocessing import LabelEncoder
  from sklearn.ensemble import RandomForestClassifier
  import scipy as stats

  #Load feature matrix and save feature names 
  df = pd.read_csv(DF, sep='\t', index_col = 0)
  #print(df[:3])
  feat_names = list(df.columns.values)[1:]

  #Recode class as 1 for positive and 0 for negative, then divide into two dataframes.
  df["Class"] = df['Class'].map({pos: 1, neg: 0})
  all_pos = df[df.Class == 1]
  pos_size = all_pos.shape[0]
  all_neg = df[df.Class == 0]
  print("%d positive examples and %d negative examples in dataframe" % (pos_size, all_neg.shape[0]))
  

  m = 0
  num_df = 3    #Number of random negative data sets to generate for balanced ML
  num_rep = 2   #Number of replicates to run with each dataset.
  num_cv = 10   #Cross validation fold

  cv_means = np.array([])
  
  #Run for each balanced dataset
  for j in range(num_df):
    
    #Make balanced dataset with random negative examples drawn
    random_neg = all_neg.sample(n=pos_size)                       #take random subset of negative examples - same size as positive
    df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe
    le=LabelEncoder()
    y = df.iloc[:,0].values                                       #Designate Column with class information
    x = df.iloc[:,1:].values                                      #Designate Columns with variable information
    cv = np.array([])
    
    #Run 'num_rep' replicates of the CV run. 
    for i in range(num_rep):                       
      m += 1
      forest = RandomForestClassifier(criterion='entropy',n_estimators=500, n_jobs=8)
      forest = forest.fit(x, y)                                                               #Train the model
      scores = cross_validation.cross_val_score(forest, x, y, cv=num_cv, scoring='f1')               #Make predictions with CV
      print(scores)
      cv = np.insert(cv, 0, np.mean(scores))
      importances = forest.feature_importances_
      if m == 1:
        imp = np.asarray([importances])             
      else:
        a = np.asarray([importances])
        imp = np.concatenate((imp, a), axis = 0)

    N = len(cv)
    F_measure_1 = np.mean(cv)
    sigma_1 = np.std(cv)
    n_1=len(cv)
    SE_1 = sigma_1/sqrt(N)
    print(cv, N, F_measure_1, sigma_1, SE_1)

    results_out.write('%s\t%.4f\t%.4f\t%.4f\n' % (j+1, F_measure_1, sigma_1, SE_1))
    print("Replicate %d: F_measure= %.4f; STDEV= %.4f; SE= %.4f" % (j+1, F_measure_1, sigma_1, SE_1))
    cv_means = np.insert(cv_means, 0, F_measure_1)

  pd_imp = pd.DataFrame(imp, index=range(1,(num_df*num_rep+1)), columns=list(full_df)[1:])   #Turn importance array into df
  pd_imp_mean = pd_imp.mean(axis=0)                    #Find means for importance measure
  pd_imp_mean.to_csv(SAVE + "_imp.txt", sep='\t')
  
  F_measure = np.mean(cv_means)
  sigma = np.std(cv_means)
  n=len(cv_means)
  SE = np.std(cv_means)/sqrt(n)
  #CI95 = stats.norm.interval(0.95, loc=F_measure, scale=sigma/np.sqrt(n))
  n_features = len(feat_names)

  results_out.write('\nNumber of features: %.1f' % (n_features))
  results_out.write('\nAverage\t%.4f\t%.4f\t%.4f' % (F_measure, sigma, SE))
  open("RESULTS.txt",'a').write('%s\t%.1f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\n' % (SAVE, n_features, num_df, num_rep, F_measure, sigma, SE))
  print('Average: F measure: %.3f\nStandard Deviation: %.3f\nStandard Error: %.3f' % (F_measure, sigma, SE))


RandomForest(DF, SAVE, pos, neg)

