"""

Merge kmer information into one database. 
Inputs:

python ~/GitHub/MotifDiscovery/pCRE_criteria.py -df ../NDD_Ras_FET01_df_p0.01.txt -enrich NDD_Ras_FET01_FETresults.txt -imp_RF NDD_Ras_FET01_RF_imp -imp_SVC NDD_Ras_FET01_SVC_imp -map ../02_kmer_mapping/all_kmers.gff

Output criteria:
pvalue - Enrichement Score 
percent_positives - Support score - i.e. Percent of cluster genes that contain the kmer in the promoter
Imp_RF - Importance score in RF predictive models
Imp_SVC - Importance score in SVC predictive models
Copy#_pos_pres - Average copy number in promoter region of cluster genes where kmer is present
Copy#_neg_pres - Average copy number in promoter region of non-cluster genes where kmer is present

Av_DHS - Proportion of kmer hits that are in accessible regions
Av_CNS - Proportion of kmer hits that are in conserved non-coding regions (between species)
Av_WSC - Average within species conservtion of kmer hits
9 - PCC_KnownTFBM - Similarity to known TFBM
10- Av_DAPSeq_Ovrp - Proportion of kmer hits that overlap with a DAP-Seq 200 bp region
11 - FS_Feature - Does the k-mer overlap with 

"""
import sys
import pandas as pd
import numpy as np
from scipy import stats


DF = ENRICH = IMP_RF = IMP_SVC = MAP = DHS = CNS = WSC = TFBM = DAP = FS = "False"

for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-df":
        DF = sys.argv[i+1]
      if sys.argv[i] == '-enrich':
        ENRICH = sys.argv[i+1]
      if sys.argv[i] == '-imp_RF':
        IMP_RF = sys.argv[i+1]
      if sys.argv[i] == "-imp_SVC":
        IMP_SVC = sys.argv[i+1]
      if sys.argv[i] == '-map':
        MAP = sys.argv[i+1]
      if sys.argv[i] == "-DHS":
        DHS = sys.argv[i+1]
      if sys.argv[i] == "-CNS":
        CNS = sys.argv[i+1]
      if sys.argv[i] == "-WSC":
        WSC = sys.argv[i+1]
      if sys.argv[i] == "-TFBM":
        TFBM = sys.argv[i+1]
      if sys.argv[i] == "-DAP_ovrp":
        DAP = sys.argv[i+1]
      if sys.argv[i] == "-FS":
        FS = int(sys.argv[i+1])


if len(sys.argv) <= 1:
  print(__doc__)
  exit()


####### Load Dataframe  #######

df_used = pd.read_csv(DF, sep='\t', header =0, index_col = 0)

kmers = np.delete(df_used.columns.values, 0)
pos_genes = df_used[df_used.Class==1].index.values
neg_genes = df_used[df_used.Class==0].index.values

df = pd.DataFrame(index = kmers, columns=None)


####### pCRE Criteria #######

if ENRICH != "False":
  enrich_df = pd.read_csv(ENRICH, sep='\t', header =0, index_col = 0, usecols = ['feature', 'percent_postives', 'pvalue'])
  df = pd.merge(df, enrich_df, how = 'left', right_index=True, left_index=True)

if IMP_RF != 'False':
  imp_RF = pd.read_csv(IMP_RF, sep='\t', header =0, index_col = 0, names = ['Imp_RF'])
  df = pd.merge(df, imp_RF, how = 'left', right_index=True, left_index=True)

if IMP_SVC != 'False':
  imp_SVC = pd.read_csv(IMP_SVC, sep='\t', header =0, index_col = 0, names = ['Imp_SVC'])
  df = pd.merge(df, imp_SVC, how = 'left', right_index=True, left_index=True)

if MAP != 'False':
  map_df = pd.read_csv(MAP, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer'])
  subset = map_df[map_df['kmer'].isin(kmers)]
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]

  count_pos = subset_pos.groupby(['kmer','gene']).size()
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).mean(), columns = ['Copy#_pos_pres'])
  num_hits_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['Hits_pos_pres'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).mean(), columns = ['Copy#_neg_pres'])
  num_hits_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['Hits_neg_pres'])

  df = pd.merge(df, grouped_pos, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, grouped_neg, how = 'left', right_index=True, left_index=True)

  hits = pd.merge(num_hits_pos, num_hits_neg, how = 'left', right_index=True, left_index=True)

if DHS != 'False':
  dhs = pd.read_csv(DHS, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  subset = dhs[dhs['kmer'].isin(kmers)]
  subset = subset[subset.details.str.contains("overlap=NA") == False]
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]


  count_pos = subset_pos.groupby(['kmer','gene']).size()
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['DHS_Hits_pos_pres'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['DHS_Hits_neg_pres'])

  hits = pd.merge(hits, grouped_pos, how = 'left', right_index=True, left_index=True)
  hits = pd.merge(hits, grouped_neg, how = 'left', right_index=True, left_index=True)
  
  hits['noDHS_Hits_pos_pres'] = hits['Hits_pos_pres'] - hits['DHS_Hits_pos_pres']
  hits['noDHS_Hits_neg_pres'] = hits['Hits_neg_pres'] - hits['DHS_Hits_neg_pres']

  hits['temp'] = hits.apply(lambda r: stats.fisher_exact([[r.DHS_Hits_pos_pres, r.noDHS_Hits_pos_pres], [r.DHS_Hits_neg_pres, r.noDHS_Hits_neg_pres]]), axis=1)
  hits = hits.replace(to_replace=['(',')'], value = ['',''], inplace = False)
  print(hits)

  #hits[['DHS_Odds', 'firstname', 'middle initial']] = df['name'].str.split(expand=True)
  #df = pd.merge(df, hits.loc[:, ['DHS(Odds,pval)']], how = 'left', right_index=True, left_index=True)

  #print(df.head(10))
