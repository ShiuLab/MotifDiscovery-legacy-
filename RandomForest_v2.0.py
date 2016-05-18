"""
PURPOSE:
Find k+mers enriched in your positive dataset and run RandomForest Classifier to determine how well those k+mers predict your classes

Before running, set path to Miniconda:  export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH

*When submitting jobs ask for 8 nodes! 

INPUT:
  -pos_file FASTA file with positive examples
  -neg_file FASTA file with negative examples - the script will run ML on 50 random samples to get balanced test
  -pos      String for what codes for the positive example (Default = 1)
  -neg      String for what codes for the negative example (Default = 0)
  -k        List of kmers to start with (/mnt/home/azodichr/ML_Python/6mers.txt or 5mers.txt)
  -pval     P-value cut off for Fisher's exact test (Default = 0.01)
  -save     Save name (will overwrite results files if the same as other names in the directory youre exporting to)
  -score    Default: F-measure. Can change to AUC-ROC using "-score roc_auc"
  -feat     Default: all (i.e. everything in the dataframe given). Can import txt file with list of features to keep.


OUTPUT:
  -SAVE_df_pPVAL.txt       Dataframe that goes into SK-learn for ML.
  -SAVE_FETresults.txt    Results for features enriched with fisher's exact test: Feature, # Pos Examples with Feature, # Neg examples with feature, Pvalue
  -SAVE_RF_results.txt    Results from RF runs
  -RESULTS.txt            Final results get added to this file: Run Name, # Features, # Reps (different Neg Datasets), CV, F_measure, StDev, SE
"""



import pandas as pd
import numpy as np
import sys
import timeit
from math import sqrt
import RF_scikit
start = timeit.default_timer()

PVAL = 0.01
FEAT = 'all'    #Features to include from dataframe. Default = all (i.e. don't remove any from the given dataframe)
SCORE = 'f1'
neg = '0'
pos = '1'

for i in range (1,len(sys.argv),2):

      if sys.argv[i] == '-pos_file':         #Fasta file for positive examples
        POS = sys.argv[i+1]
      if sys.argv[i] == '-neg_file':         #Fasta file for negative examples
        NEG = sys.argv[i+1]
      if sys.argv[i] == '-neg':              #String for negative class : Default = 0
        neg = sys.argv[i+1]
      if sys.argv[i] == "-pos":              #String for positive class : Default = 1
        pos = sys.argv[i+1]
      if sys.argv[i] == '-pval':             #Default is 0.01
        PVAL = float(sys.argv[i+1])
      if sys.argv[i] == '-save':
        SAVE = sys.argv[i+1]
      if sys.argv[i] == '-feat':
        FEAT = sys.argv[i+1]
      if sys.argv[i] == '-k':
        K = sys.argv[i+1]
      if sys.argv[i] == "-score":
        SCORE = sys.argv[i+1]


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

  DF= pd.DataFrame(dataframe, index=genes, columns=numpy_header, dtype=int)  # , dtype=int # Turn numpy into pandas DF
  DF= DF.drop("Skip_this_line",0)

  
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
        DF = DF.drop(k, 1)
      if count%10000==0:
        print("Completed " + str(count) + " features")

    except ValueError:
      count += 1 
      outFISH.write('\n%s\t%d\t%d\t1.0' % (k, (positive_present[k]),(negative_present[k])))
  return(DF, enriched_kmers)



def Make_DF(K, PVAL, SAVE):
  print("Testing all possible kmers given")
  #Put all kmers/kmer pairs into list
  kmer_5 = []
  for l in open(K, 'r'):
    kmer_5.append(l.strip("\n"))

  DF_temp, enriched_5mers = Find_Enrich(POS, NEG, kmer_5, PVAL, SAVE)
  final_km = []

  print("Testing %d k+1mers" % (len(enriched_5mers)*8))
  kmer_6 = {}
  for key in enriched_5mers:
    final_km.append(key)
    temp_km = []
    temp_km.extend([key+"A", key+"T", key+"G", key+"C", "A"+key, "T"+key, "G"+key, "C"+key])
    DF_temp, enriched_6mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
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
    DF_temp, enriched_7mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
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
    DF_temp, enriched_8mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
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
    DF_temp, enriched_9mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
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
    DF_temp, enriched_10mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
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
    DF_temp, enriched_11mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
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
    DF_temp, enriched_12mer = Find_Enrich(POS, NEG, temp_km, PVAL, SAVE)
    if len(enriched_12mer) > 0:   
      for k in enriched_12mer:
        if enriched_12mer[k] <= kmer_11[key]:
          kmer_12[k]=enriched_12mer[k]        #If 12mer has lower pvalue than the 11mer, add it to a list to retest as a 7mer.
          final_km.append(k)    

  final_km_NoDups = list(set(final_km))
  DF, enriched_kmers = Find_Enrich(POS, NEG, final_km_NoDups, PVAL, SAVE)
  n_features = DF.shape[1]-2
  print('%d kmers enriched at p = %s' % (n_features, PVAL))
  enriched_name = SAVE + "_df_p" + str(PVAL) + ".txt"
  DF["Class"] = DF["Class"].replace(1, pos)
  DF["Class"] = DF["Class"].replace(0, neg)
  DF.to_csv(enriched_name, sep='\t')
  return(DF)


if __name__ == '__main__':
  
  RF_scikit.RandomForest(Make_DF(K, PVAL, SAVE), SAVE, SCORE, FEAT, pos, neg)

stop = timeit.default_timer()
print('Run time: %.2f min' % ((stop-start)/60))

