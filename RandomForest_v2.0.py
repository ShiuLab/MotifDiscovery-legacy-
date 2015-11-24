import pandas as pd
import numpy as np
import sys
import timeit
from math import sqrt
start = timeit.default_timer()

IMP = 'no'  #Must specifiy -importance yes if you want importance calculated
PVAL = 0.05

for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-df":
        pd_DF = sys.argv[i+1]
      if sys.argv[i] == '-imp':         #Set to 'yes' if you want importance values calculated
        IMP = sys.argv[i+1]
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
  
  dataframe = np.zeros([1,len(km)+1])       # This is a temporary fix for fitting the np df into the pd df - the plus 1 is for the Class!
  positive_present = {}.fromkeys(km, 0)     # Count occurence of each feature in positive examples
  negative_present = {}.fromkeys(km, 0)     # Count occurence of each feature in negative examples

  #Open positive and negative fasta files
  p = open(POS, 'r')
  n = open(NEG, 'r')

  #print("Making dataframe...")
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
        if str(k1) in seq or str(k1.reverse_complement()) in seq and str(k2) in seq or str(k2.reverse_complement()) in seq:  #If k1 and k2 in seq
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
  #print("Positive examples in dataframe")
  
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
  #print("Negative examples in dataframe")

  pd_DF= pd.DataFrame(dataframe, index=genes, columns=numpy_header, dtype=int)  # Turn numpy into pandas DF
  pd_DF= pd_DF.drop("Skip_this_line",0)

  
  # Calculate enrichement scores
  #print("Running Fisher's Exact test...")
  outFISH = open(SAVE+"_FETresults.txt",'w')
  outFISH.write('feature\tPosCount\tNegCount\tpvalue')
  count = 0
  
  enriched_kmers = {}
  for k in positive_present:
    try:
      count += 1 
      oddsratio,pvalue = fisher_exact([[positive_present[k],(num_pos-positive_present[k])],[negative_present[k],(num_neg-negative_present[k])]])
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
  print("Testing all possible 5mers")
  #Put all kmers/kmer pairs into list
  kmer_5 = []
  for l in open(K, 'r'):
    kmer_5.append(l.strip("\n"))

  pd_DF_temp, enriched_5mers = Find_Enrich(POS, NEG, kmer_5, PVAL, SAVE)
  final_km = []

  print("Testing %d 6mers" % (len(enriched_5mers)*8))
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

  print("Testing %d 7mers" % (len(kmer_6)*8))
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

  print("Testing %d 8mers" % (len(kmer_7)*8))
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

  print("Testing %d 9mers" % (len(kmer_8)*8))
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
  
  print("Testing %d 10mers" % (len(kmer_9)*8))
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

  print("Testing %d 11mers" % (len(kmer_10)*8))
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

  print("Testing %d 12mers" % (len(kmer_11)*8))
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

  print('%d kmers enriched at p = %s' % (pd_DF.shape[1]-1, PVAL))
  enriched_name = SAVE + "_df_p" + str(PVAL) + ".txt"
  pd_DF.to_csv(enriched_name, sep='\t')
  return(pd_DF)


def RandomForest(pd_DF, SAVE):

  from sklearn.cross_validation import train_test_split, cross_val_score, ShuffleSplit
  from sklearn.ensemble import RandomForestClassifier 

  #full_df = pd.read_csv(pd_DF, sep='\t',header=0)
  #df = pd_DF
  #var_names = list(df.columns.values)[2:]

  full_df = pd_DF
  var_names = list(full_df.columns.values)[2:]
  results_out = open(SAVE + "_RF_results.txt",'w')
  
  all_pos = full_df[full_df.Class == 1]
  pos_size = all_pos.shape[0]
  all_neg = full_df[full_df.Class == 0]
  print("%d positive examples and %d negative examples in dataframe" % (pos_size, all_neg.shape[0]))
  
  m = 0

  num_df_to_run = 20
  num_cv_to_run = 10

  cv_means = np.array([])
  
  for j in range(num_df_to_run):
    
    random_neg = all_neg.sample(n=pos_size)               #take random subset of negative examples - same size as positive
    df = pd.DataFrame.merge(all_pos, random_neg, how = 'outer')  #Make balanced dataframe
    y = df.iloc[:,0].values                               #Designate Column with class information
    x = df.iloc[:,1:].values                              #Designate Columns with variable information
    cv = np.array([])
    
    for i in range(num_cv_to_run):                        #Run x replicates of cv
      m += 1
      X_train, X_test, Y_train, Y_test = train_test_split(x,y,test_size=0.3, random_state=i)              #Set up training/testing data sets for each replication
      forest = RandomForestClassifier(criterion='entropy',n_estimators=500, random_state=1, n_jobs=2)
      forest = forest.fit(X_train, Y_train)                                                               #Train the model
      scores = cross_val_score(estimator=forest, X=X_train, y=Y_train, cv=10, n_jobs=2)                   #Make predictions with 10x CV
      cv = np.insert(cv, 0, np.mean(scores))
      if IMP == 'yes':
        importances = forest.feature_importances_
        if m == 1:
          imp = np.asarray([importances])             
        else:
          a = np.asarray([importances])
          imp = np.concatenate((imp, a), axis = 0)
    results_out.write('Random Dataframe %s: CV accuracy: %.3f 95%% CI +/- %.3f\n' % (j+1, np.mean(cv), np.std(cv)*2))
    print("Replicate %d: CV accuracy: %.3f 95%% CI +/- %.3f" % (j+1, np.mean(cv), np.std(cv)*2))
    cv_means = np.insert(cv_means, 0, np.mean(cv))

  if IMP == 'yes':
    pd_imp = pd.DataFrame(imp, index=range(1,(num_df_to_run*num_cv_to_run+1)), columns=list(full_df)[1:])   #Turn importance array into df
    #pd_imp.to_csv(SAVE + "_imp_detailed.txt", sep='\t', header = "Feature\tImportanceMeasure")              #Save details - can turn this off!
    pd_imp_mean = pd_imp.mean(axis=0)                    #Find means for importance measure
    pd_imp_mean.to_csv(SAVE + "_imp.txt", sep='\t')
    
  results_out.write('\nAverage: CV accuracy: %.3f 95%% CI +/- %.3f\n' % (np.mean(cv_means), np.std(cv_means)*2))
  open("RESULTS.txt",'a').write('%s\t%.3f\t%.3f\t%.3f\n' % (SAVE, np.mean(cv_means), np.std(cv_means)*2,(np.std(cv_means)/sqrt(num_df_to_run))))
  print('Average: CV accuracy: %.3f 95%% CI +/- %.3f\nStandard Deviation: %.3f\nStandard Error: %.3f' % (np.mean(cv_means), np.std(cv_means)*2, np.std(cv_means), (np.std(cv_means)/sqrt(num_df_to_run))))

#Make_DF(K, PVAL, SAVE)
#RandomForest(pd_DF,SAVE)
#RandomForest(Make_DF(POS, NEG, K, PVAL, SAVE),SAVE)

RandomForest(Make_DF(K, PVAL, SAVE), SAVE)

stop = timeit.default_timer()
print('Run time: %.2f min' % ((stop-start)/60))

"""
#Other options available for ML using sklearn


# Get Confusion Matrix from one run
from sklearn.metrics import confusion_matrix
y_pred = forest.predict(X_test)
confmat = confusion_matrix(y_true=Y_test, y_pred=y_pred)


#Just one rep with no CV - print Precision, Recall and F score
from sklearn.metrics import precision_score, recall_score, f1_score
print('Precision: %.3f' % precision_score(y_true=Y_test, y_pred=y_pred))
print('Recall: %.3f' % recall_score(y_true=Y_test, y_pred=y_pred))
print('F1: %.3f' % f1_score(y_true=Y_test, y_pred=y_pred))

#Long way to calculate CV score
from sklearn.cross_validation import StratifiedKFold
kfold = StratifiedKFold(y=Y_train, n_folds=20, random_state=1)
scores = []
for k, (train, test) in enumerate(kfold):
  forest.fit(X_train[train],Y_train[train])
  score = forest.score(X_train[test],Y_train[test])
  scores.append(score)
  print('Fold: %s, Class dist: %s, Acc: %.3f' % (k+1, np.bincount(Y_train[train]),score))
print('CV accuracy: %.3f +/- %.3f' % (np.mean(scores), np.std(scores)))


# To plot confusion matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax=plt.subplots(2,2,figsize=(2.5,2.5))
ax.matshow(confmat, cmap=plt.cm.Blues, alpha=0.3)
for i in range(confmat.shape[0]):
  for j in range(confmat.shape[1]):
    ax.text(x=j, y=i, s=confmat[i,j], va='center', ha='center')
plt.xlabel('predicted')
plt.ylabel('true')
plt.savefig("test_fig.png")


#To plot decision if numerical variables (i.e. petal width vs. petal length)
from plot_decision import plot_decision_regions
X_combined = np.vstack((X_train, X_test))
Y_combined = np.hstack((Y_train, Y_test))
plot_decision_regions(X_combined, Y_combined, classifier=forest, test_idx=range(105,150))
"""
