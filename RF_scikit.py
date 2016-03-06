"""
PURPOSE:
Run classification RF from sci-kit learn on a given dataframe
*When submitting jobs ask for 8 nodes! 

Must set path to Miniconda in HPC:  export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH


INPUT:
  -df       Feature dataframe for ML. Format -> Col 1 = example.name, Col 2 = Class, Col 3-... = Features.
  -save     Save name (will overwrite some results files if the same as other names in the directory youre exporting to)
  -pos      String for what codes for the positive example (i.e. UUN) Default = 1
  -neg      String for what codes for the negative example (i.e. NNN) Default = 0
  -score    Default: F-measure. Can change to AUC-ROC using roc_auc
  -feat     Default: all (i.e. everything in the dataframe given). Can import txt file with list of features to keep.

OUTPUT:
  -SAVE_RF_results.txt    Results from RF runs
  -SAVE_imp.txt           Importance scores for each feature
  -RESULTS.txt            Final results get added to this file: Run Name, #Pos Examples, #Neg Examples, # Features, # Reps (different Neg Datasets), #CV Reps, CV Fold, F_measure, StDev, SE
"""

import pandas as pd
import numpy as np
import sys
from math import sqrt

neg = int(0)
pos = int(1)
SCORE = 'f1'    #Scoring method for RF, default F-measure, can change to AUC-ROC using -score roc_auc
FEAT = 'all'    #Features to include from dataframe. Default = all (i.e. don't remove any from the given dataframe)

for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-score":
        SCORE = sys.argv[i+1]
      if sys.argv[i] == "-df":
        DF = sys.argv[i+1]
      if sys.argv[i] == '-save':
        SAVE = sys.argv[i+1]
      if sys.argv[i] == '-feat':
        FEAT = sys.argv[i+1]
      if sys.argv[i] == '-neg':
        neg = sys.argv[i+1]
      if sys.argv[i] == "-pos":
        pos = sys.argv[i+1]


def RandomForest(DF, SAVE, SCORE, FEAT, pos, neg):

  from sklearn import cross_validation
  from sklearn.ensemble import RandomForestClassifier
  import scipy as stats

  print(type(DF))
  #Load feature matrix and save feature names 
  if isinstance(DF, str):
    df = pd.read_csv(DF, sep='\t', index_col = 0)
  else:
    df = DF

  #If features to keep list given, remove columns not in that list
  if FEAT != 'all':
    with open(FEAT) as f:
      features = f.read().splitlines()
      features = ['Class'] + features
    df = df.loc[:,features]
  
  feat_names = list(df.columns.values)[1:]
  print(df)
  #Recode class as 1 for positive and 0 for negative, then divide into two dataframes.
  df["Class"] = df['Class'].map({pos: 1, neg: 0})
  print(df)
  all_pos = df[df.Class == 1]
  pos_size = all_pos.shape[0]
  all_neg = df[df.Class == 0]
  print("%d positive examples and %d negative examples in dataframe" % (pos_size, all_neg.shape[0]))
  
  #Set output file
  results_out = open(SAVE + "_RF_results.txt",'w')
  results_out.write('Dataframe_Rep\tF_measure\tSTDEV\tSE')

  m = 0 
  num_df = 50        #Number of random balanced replicates (Rand_df_reps)
  num_rep = 10       #Number of CV replicates (cv_reps)
  num_cv = 10       #Cross validation fold


  #Make empty array to save score from each random balanced replicate. Size = num_df
  Rand_df_reps = np.array([])
  
  #Run for each balanced dataset
  for j in range(num_df):
    
    #Make balanced dataset with random negative examples drawn
    random_neg = all_neg.sample(n=pos_size)                       #take random subset of negative examples - same size as positive
    df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')   #Make balanced dataframe
    y = df.iloc[:,0].values                                       #Designate Column with class information
    x = df.iloc[:,1:].values                                      #Designate Columns with variable information
    
    #Make an empty array to save scores each CV replicate. Size = num_rep
    cv_reps = np.array([])
    
    #Run 'num_rep' replicates of the CV run. 
    for i in range(num_rep):                       
      m += 1
      forest = RandomForestClassifier(criterion='entropy',n_estimators=500, n_jobs=8)         #Define the type of model
      forest = forest.fit(x, y)                                                               #Train the model
      scores = cross_validation.cross_val_score(forest, x, y, cv=num_cv, scoring=SCORE)     #Make predictions with CV
      cv_reps = np.insert(cv_reps, 0, np.mean(scores))

      #Add importance values to imp array
      importances = forest.feature_importances_
      if m == 1:
        imp = np.asarray([importances])             
      else:
        a = np.asarray([importances])
        imp = np.concatenate((imp, a), axis = 0)

    #Calculate stats for the CV replicates, write to _RF_results.txt file and print
    F_measure_cv_reps = np.mean(cv_reps)
    sigma_cv_reps = np.std(cv_reps)
    SE_cv_reps = sigma_cv_reps/sqrt(len(cv_reps))
    results_out.write('\n%s\t%.4f\t%.4f\t%.4f' % (j+1, F_measure_cv_reps, sigma_cv_reps, SE_cv_reps))
    print("Replicate %d: F_measure= %.4f; STDEV= %.4f; SE= %.4f" % (j+1, F_measure_cv_reps, sigma_cv_reps, SE_cv_reps))
    
    #Add score mean from the random balanced dataset replicate to the Rand_df_reps array. 
    Rand_df_reps = np.insert(Rand_df_reps, 0, F_measure_cv_reps)

  #Output mean of the importance measure for each feature
  pd_imp = pd.DataFrame(imp, index=range(1,(num_df*num_rep+1)), columns=feat_names)   #Turn importance array into df
  pd_imp_mean = pd_imp.mean(axis=0)                    #Find means for importance measure
  pd_imp_mean.to_csv(SAVE + "_imp.txt", sep='\t')
  
  #Calculate stats for the random dataframe replicates, write to the bottom of _RF_results.txt, add to RESULTS.txt, and print
  F_measure = np.mean(Rand_df_reps)
  sigma = np.std(Rand_df_reps)
  SE = np.std(Rand_df_reps)/sqrt(num_df)
  #CI95 = stats.norm.interval(0.95, loc=F_measure, scale=sigma/np.sqrt(n))
  n_features = len(feat_names)

  results_out.write('\nNumber of features: %.1f' % (n_features))
  results_out.write('\nAverage\t%.4f\t%.4f\t%.4f' % (F_measure, sigma, SE))
  open("RESULTS.txt",'a').write('%s\t%i\t%i\t%i\t%i\t%i\t%i\t%.5f\t%.5f\t%.5f\n' % (SAVE, pos_size, all_neg.shape[0], n_features, num_df, num_rep, num_cv, F_measure, sigma, SE))
  print('\nAverage: F measure: %.3f; stdev: %.3f; SE: %.3f' % (F_measure, sigma, SE))
  print ('\nColumn Names in RESULTS.txt output: Run_Name, #Pos_Examples, #Neg_Examples, #Features, #Random_df_Reps, #CV_Reps, CV_Fold, F_measure, StDev, SE')


if __name__ == "__main__":
  RandomForest(DF, SAVE, SCORE, FEAT, pos, neg)

