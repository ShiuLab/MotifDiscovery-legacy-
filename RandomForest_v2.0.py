"""
PURPOSE:
Find k+mers enriched in your positive dataset and run RandomForest Classifier to determine how well those k+mers predict your classes

*When submitting jobs ask for 8 nodes! 

INPUT:
  -pos      FASTA file with positive examples
  -neg      FASTA file with negative examples - can be much larger, the script will randomly sample for each ML run to get balanced test
  -k        List of kmers to start with (/mnt/home/azodichr/ML_Python/6mers.txt or 5mers.txt)
  -pval     P-value cut off for Fisher's exact test (Default = 0.05)
  -save     Save name (will overwrite some results files if the same as other names in the directory youre exporting to)

OUTPUT:
  -SAVE_df_pXXX.txt       Dataframe that goes into SK-learn for ML.
  -SAVE_FETresults.txt    Results for features enriched with fisher's exact test: Feature, # Pos Examples with Feature, # Neg examples with feature, Pvalue
  -SAVE_RF_results.txt    Results from RF runs
  -RESULTS.txt            Final results get added to this file: Run Name, # Features, # Reps (different Neg Datasets), CV, F_measure, StDev, SE
"""



import pandas as pd
import numpy as np
import sys
import timeit
from math import sqrt
start = timeit.default_timer()

PVAL = 0.05

for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-df":
        pd_DF = sys.argv[i+1]
      if sys.argv[i] == '-pos':         #Fasta file for positive examples
        POS = sys.argv[i+1]
      if sys.argv[i] == '-neg':         #Fasta file for negative examples
        NEG = sys.argv[i+1]
      if sys.argv[i] == '-pval':        #Default is 0.05
        PVAL = float(sys.argv[i+1])
      if sys.argv[i] == '-save':
        SAVE = sys.argv[i+1]
      if sys.argv[i] == '-k':
        K = sys.argv[i+1]


def Find_Enrich(POS, NEG, km, PVAL, SAVE):
  from Bio import SeqIO
  from Bio.Seq import Seq
  from scipy.stats import fisher_exact
  
  numpy_header = ['Class']
  
  for i in km:
    numpy_header.append(i)  
  
  dataframe = np.zeros([1,len(km)+1])       # This fits the np df into the pd df - the plus 1 is for the Class!
  positive_present = {}.fromkeys(km, 0)     # Count occurence of each feature in positive examples
  negative_present = {}.fromkeys(km, 0)     # Count occurence of each feature in negative examples

  #Open positive and negative fasta files
  p = open(POS, 'r')
  n = open(NEG, 'r')

  num_pos = 0
  num_neg = 0
  genes = ['Skip_this_line']  #index for pandas df

  for seq_record in SeqIO.parse(p, 'fasta'):
    num_pos += 1
    header = seq_record.id
    genes.append(header)
    seq = str(seq_record.seq)
    gene_array =np.array([1])       # Array of P/A (1/0) for each gene - starts with '1' For Positive Class
    for ki in km:
      if " " in ki:                         #Checks to see if motif is a pair - pairs are separated by a space
        k1 = Seq(ki.split(" ")[0])
        k2 = Seq(ki.split(" ")[1])
        if str(k1) in seq or str(k1.reverse_complement()) in seq and str(k2) in seq or str(k2.reverse_complement()) in seq: 
          gene_array = np.append(gene_array, 1)
          positive_present[ki] = positive_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)

      else:                                #If no separation by a space, assumes you're looking at singletons.
        kmer = Seq(ki)
        if str(kmer) in seq or str(kmer.reverse_complement()) in seq:
          gene_array = np.append(gene_array, 1)
          positive_present[ki] = positive_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)
    dataframe = np.vstack((dataframe,gene_array))
  
  for seq_record in SeqIO.parse(n, 'fasta'):
    num_neg += 1 
    header = seq_record.id
    genes.append(header)
    seq = str(seq_record.seq)         
    gene_array =np.array([0])       # Array of P/A (1/0) for each gene - starts with '0' For Negative Class
    for ki in km:
      if " " in ki:                         #Checks to see if motif is a pair - pairs are separated by a space
        k1 = Seq(ki.split(" ")[0])
        k2 = Seq(ki.split(" ")[1])
        if str(k1) in seq or str(k1.reverse_complement()) in seq and str(k2) in seq or str(k2.reverse_complement()) in seq:
          gene_array = np.append(gene_array, 1)
          negative_present[ki] = negative_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)

      else:                                 #If no separation by a space, assumes you're looking at singletons.
        kmer = Seq(ki)
        if str(kmer) in seq or str(kmer.reverse_complement()) in seq:
          gene_array = np.append(gene_array, 1)
          negative_present[ki] = negative_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)
    dataframe = np.vstack((dataframe,gene_array))

  pd_DF= pd.DataFrame(dataframe, index=genes, columns=numpy_header, dtype=int)  # Turn numpy into pandas DF
  pd_DF= pd_DF.drop("Skip_this_line",0)

  
  # Calculate enrichement scores
  outFISH = open(SAVE+"_FETresults.txt",'w')
  outFISH.write('feature\tPosCount\tNegCount\tpvalue')
  count = 0
  
  enriched_kmers = {}
  for k in positive_present:
    try:
      count += 1
      TP = positive_present[k]            #Positive Examples with kmer present
      FP = negative_present[k]            #Negative Examples with kmer present
      TN = num_neg-negative_present[k]    #Negative Examples without kmer
      FN = num_pos-positive_present[k]    #Positive Examples without kmer

      oddsratio,pvalue = fisher_exact([[TP,FN],[FP,TN]],alternative='greater')
      outFISH.write('\n%s\t%d\t%d\t%.7f' % (k, (positive_present[k]),(negative_present[k]),pvalue))
      if pvalue <= PVAL:          # Remove unenriched features from dataframe
        enriched_kmers[k] = pvalue
      if pvalue > PVAL:
        pd_DF = pd_DF.drop(k, 1)
      if count%10000==0:
        print("Completed " + str(count) + " features")

    except ValueError:
      count += 1 
      outFISH.write('\n%s\t%d\t%d\t1.0' % (k, (positive_present[k]),(negative_present[k])))
  return(pd_DF, enriched_kmers)



def Make_DF(K, PVAL, SAVE):
  print("Testing all possible kmers given")
  #Put all kmers/kmer pairs into list
  kmer_5 = []
  for l in open(K, 'r'):
    kmer_5.append(l.strip("\n"))

  pd_DF_temp, enriched_5mers = Find_Enrich(POS, NEG, kmer_5, PVAL, SAVE)
  final_km = []

  print("Testing %d k+1mers" % (len(enriched_5mers)*8))
  kmer_6 = {}
  for key in enriched_5mers:
    final_km.append(key)
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    pd_DF_temp, enriched_6mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_6mer) > 0:               #If enriched see if the pval is lower than it was for the 5mer
      for k in enriched_6mer:
        if enriched_6mer[k] <= enriched_5mers[key]:
          kmer_6[k]=enriched_6mer[k]          #If 6mer has lower pvalue than the 5mer, add it to a list to retest as a 7mer.
          final_km.append(k) 

  print("Testing %d k+2mers" % (len(kmer_6)*8))
  kmer_7 = {}
  for key in kmer_6:
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    pd_DF_temp, enriched_7mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_7mer) > 0:   
      for k in enriched_7mer:
        if enriched_7mer[k] <= kmer_6[key]:
          kmer_7[k]=enriched_7mer[k]        #If 7mer has lower pvalue than the 6mer, add it to a list to retest as a 8mer.
          final_km.append(k)

  print("Testing %d k+3mers" % (len(kmer_7)*8))
  kmer_8 = {}
  for key in kmer_7:
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    pd_DF_temp, enriched_8mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_8mer) > 0:   
      for k in enriched_8mer:
        if enriched_8mer[k] <= kmer_7[key]:
          kmer_8[k]=enriched_8mer[k]        #If 8mer has lower pvalue than the 7mer, add it to a list to retest as a 9mer.
          final_km.append(k)

  print("Testing %d k+4mers" % (len(kmer_8)*8))
  kmer_9 = {}
  for key in kmer_8:
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    pd_DF_temp, enriched_9mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_9mer) > 0:   
      for k in enriched_9mer:
        if enriched_9mer[k] <= kmer_8[key]:
          kmer_9[k]=enriched_9mer[k]        #If 9mer has lower pvalue than the 8mer, add it to a list to retest as a 10mer.
          final_km.append(k)
  
  print("Testing %d k+5mers" % (len(kmer_9)*8))
  kmer_10 = {}
  for key in kmer_9:
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    pd_DF_temp, enriched_10mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_10mer) > 0:   
      for k in enriched_10mer:
        if enriched_10mer[k] <= kmer_9[key]:
          kmer_10[k]=enriched_10mer[k]        #If 10mer has lower pvalue than the 9mer, add it to a list to retest as a 11mer.
          final_km.append(k)

  print("Testing %d k+6mers" % (len(kmer_10)*8))
  kmer_11 = {}
  for key in kmer_10:
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    pd_DF_temp, enriched_11mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_11mer) > 0:   
      for k in enriched_11mer:
        if enriched_11mer[k] <= kmer_10[key]:
          kmer_11[k]=enriched_11mer[k]        #If 11mer has lower pvalue than the 10mer, add it to a list to retest as a 12mer.
          final_km.append(k)  

  print("Testing %d k+7mers" % (len(kmer_11)*8))
  kmer_12 = {}
  for key in kmer_11:
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    pd_DF_temp, enriched_12mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_12mer) > 0:   
      for k in enriched_12mer:
        if enriched_12mer[k] <= kmer_11[key]:
          kmer_12[k]=enriched_12mer[k]        #If 12mer has lower pvalue than the 11mer, add it to a list to retest as a 7mer.
          final_km.append(k)    

  final_km_NoDups = list(set(final_km))
  pd_DF, enriched_kmers = Find_Enrich(POS, NEG, final_km_NoDups, PVAL, SAVE)
  n_features = pd_DF.shape[1]-2
  print('%d kmers enriched at p = %s' % (n_features, PVAL))
  enriched_name = SAVE + "_df_p" + str(PVAL) + ".txt"
  pd_DF.to_csv(enriched_name, sep='\t')
  return(pd_DF)


def RandomForest(pd_DF, SAVE):

  from sklearn.cross_validation import train_test_split, cross_val_score, ShuffleSplit
  from sklearn.ensemble import RandomForestClassifier 
  import scipy as stats

  #full_df = pd.read_csv(pd_DF, sep='\t',header=0)
  #df = pd_DF
  #var_names = list(df.columns.values)[2:]

  full_df = pd_DF
  var_names = list(full_df.columns.values)[2:]
  results_out = open(SAVE + "_RF_results.txt",'w')
  results_out.write('Dataframe_Rep\tF_measure\tSTDEV\tSE\n')

  all_pos = full_df[full_df.Class == 1]
  pos_size = all_pos.shape[0]
  all_neg = full_df[full_df.Class == 0]
  print("%d positive examples and %d negative examples in dataframe" % (pos_size, all_neg.shape[0]))
  
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

#Make_DF(K, PVAL, SAVE)
#RandomForest(pd_DF,SAVE)

RandomForest(Make_DF(K, PVAL, SAVE), SAVE)

stop = timeit.default_timer()
print('Run time: %.2f min' % ((stop-start)/60))

