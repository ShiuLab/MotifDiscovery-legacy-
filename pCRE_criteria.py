"""

Merge kmer information into one database. 
Inputs:

python ~/GitHub/MotifDiscovery/pCRE_criteria.py -df ../NDD_Ras_FET01_df_p0.01.txt -enrich NDD_Ras_FET01_FETresults.txt -imp_RF NDD_Ras_FET01_RF_imp -imp_SVC NDD_Ras_FET01_SVC_imp -map ../02_kmer_mapping/all_kmers.gff -DHS ../02_kmer_mapping/all_kmers.gff.DHS_all_ovrlp -CNS ../02_kmer_mapping/all_kmers.gff.CNS_between_ovrlp  -FS NDD_Ras_df.txt_decisiontree_50 -TFBM_DAPSeq ../03_KnownTFBM_Similarity/02_DAP_Seq/kmers.txt.tamo-DAP_motifs.txt.tm_mod.txt -TFBM_CISBP ../03_KnownTFBM_Similarity/01_CIS_BP/kmers.txt.tamo-Athaliana_TFBM_v1.01.tm.index.direct.index.tm_mod.txt

##### Input #####
Required:
-df                 Dataframe with Col 1 = genes, Col 2 = Class, Col 3-... = pCREs

Optional:
-enrich             FET results from pCRE finding pipeline (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/01_Kmer_Features/NDD_Ras_FET01_FETresults.txt)
-imp_RF             RF Machine learning importance file 
-imp_SVC            SVC Machine learning importance file
-map                GFF file of all pCREs of interest (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/02_kmer_mapping/all_kmers.gff)
-DHS                Also required (-map)
-DAP                
-CNS                
-WSC                
-TFBM_CISBP         PCCD between kmers and CIS-BP motifs (~mjliu/kmer_5/Athaliana_TFBM_v1.01.tm.index.direct.index.tm) using output from Ming's pipeline (~mjliu/kmer_5/kmer_5.sh)
-TFBM_DAPSeq        PCCD between kmers and DAP-Seq motifs (/mnt/research/ShiuLab/14_DAPseq/PWM_to_tamo/DAP_motifs.txt.tm) using output from Ming's pipeline (~mjliu/kmer_5/kmer_5.sh)
-FS                 pCREs selected by decision tree feature selection rather than enrichement.

##### Possible Output Columns #####
pvalue              Enrichement Score 
percent_positives   Support score - i.e. Percent of cluster genes that contain the kmer in the promoter
Imp_RF              Importance score in RF predictive models
Imp_SVC             Importance score in SVC predictive models
Copy#_pos_pres      Average copy number in promoter region of cluster genes where kmer is present
Copy#_neg_pres      Average copy number in promoter region of non-cluster genes where kmer is present
DHS_Odds            Fisher's Exact Test comparing pCREs overlapping with DHS sites in cluster and non-cluster genes
DHS_pval            Significance score for DHS FET
CNS_Odds            Fisher's Exact Test comparing pCREs overlapping with CNS sites in cluster and non-cluster genes
CNS_pval            Significance score for CNS FET
FS_ovlp             Does the k-mer overlap with 
CISBP_TFBM          Best match known TFBM from the CIS-BP Database
CISBP_TFBM_PCCD     Similarity (PCC-Distance - lower = better match) to best match known TFBM from the CIS-BP Database
DAPSeq_TFBM         Best match known TFBM from DAP-Seq Database
DAPSeq_TFBM_PCCD    Similarity (PCC-Distance - lower = better match) to best match known TFBM from DAP-Seq Database


Av_WSC - Average within species conservtion of kmer hits
DAPSeq_Odds - Fisher's Exact Test comparing pCREs overlapping with DAP-Seq hit sites (200 bp regions) in cluster and non-cluster genes
DAPSeq_pval - Significance score for DAPSeq FET


Fisher's Exact Test
                             Overlap with Feature   |  No overlap with feature
pCRE in Cluster gene      |       t1                            t2
pCRE in Non-Cluster gene  |       t3                            t4

"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

DF = ENRICH = IMP_RF = IMP_SVC = MAP = DHS = CNS = WSC = TFBM = DAP = FS = TFBM_CISBP = TFBM_DAPSeq = False

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
      if sys.argv[i] == "-TFBM_CISBP":
        TFBM_CISBP = sys.argv[i+1]
      if sys.argv[i] == "-TFBM_DAPSeq":
        TFBM_DAPSeq = sys.argv[i+1]
      if sys.argv[i] == "-DAP":
        DAP = sys.argv[i+1]
      if sys.argv[i] == "-FS":
        FS = sys.argv[i+1]


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

if ENRICH:
  enrich_df = pd.read_csv(ENRICH, sep='\t', header =0, index_col = 0, usecols = ['feature', 'percent_postives', 'pvalue'])
  df = pd.merge(df, enrich_df, how = 'left', right_index=True, left_index=True)

if IMP_RF:
  imp_RF = pd.read_csv(IMP_RF, sep='\t', header =0, index_col = 0, names = ['Imp_RF'])
  df = pd.merge(df, imp_RF, how = 'left', right_index=True, left_index=True)

if IMP_SVC:
  imp_SVC = pd.read_csv(IMP_SVC, sep='\t', header =0, index_col = 0, names = ['Imp_SVC'])
  df = pd.merge(df, imp_SVC, how = 'left', right_index=True, left_index=True)

if MAP:
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

if DHS:
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
  hits = hits.fillna(0)
  hits['temp'] = hits.apply(lambda r: stats.fisher_exact([[r.DHS_Hits_pos_pres, r.noDHS_Hits_pos_pres], [r.DHS_Hits_neg_pres, r.noDHS_Hits_neg_pres]]), axis=1)
  hits[['DHS_Odds', 'DHS_pval']] = hits['temp'].apply(pd.Series)

  df = pd.merge(df, hits.loc[:, ['DHS_Odds', 'DHS_pval']], how = 'left', right_index=True, left_index=True)

if CNS:
  cns = pd.read_csv(CNS, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  subset = cns[cns['kmer'].isin(kmers)]
  subset = subset[subset.details.str.contains("overlap=NA") == False]
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]


  count_pos = subset_pos.groupby(['kmer','gene']).size()
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['CNS_Hits_pos_pres'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['CNS_Hits_neg_pres'])

  hits = pd.merge(hits, grouped_pos, how = 'left', right_index=True, left_index=True)
  hits = pd.merge(hits, grouped_neg, how = 'left', right_index=True, left_index=True)
  
  hits['noCNS_Hits_pos_pres'] = hits['Hits_pos_pres'] - hits['CNS_Hits_pos_pres']
  hits['noCNS_Hits_neg_pres'] = hits['Hits_neg_pres'] - hits['CNS_Hits_neg_pres']
  hits = hits.fillna(0)
  hits['temp'] = hits.apply(lambda r: stats.fisher_exact([[r.CNS_Hits_pos_pres, r.noCNS_Hits_pos_pres], [r.CNS_Hits_neg_pres, r.noCNS_Hits_neg_pres]]), axis=1)
  hits[['CNS_Odds', 'CNS_pval']] = hits['temp'].apply(pd.Series)
  
  df = pd.merge(df, hits.loc[:, ['CNS_Odds', 'CNS_pval']], how = 'left', right_index=True, left_index=True)

if DAP:
  dap = pd.read_csv(DAP, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  subset = dap[dap['kmer'].isin(kmers)]
  subset = subset[subset.details.str.contains("overlap=NA") == False]
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]


  count_pos = subset_pos.groupby(['kmer','gene']).size()
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['DAP_Hits_pos_pres'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['DAP_Hits_neg_pres'])

  hits = pd.merge(hits, grouped_pos, how = 'left', right_index=True, left_index=True)
  hits = pd.merge(hits, grouped_neg, how = 'left', right_index=True, left_index=True)
  
  hits['noDAP_Hits_pos_pres'] = hits['Hits_pos_pres'] - hits['DAP_Hits_pos_pres']
  hits['noDAP_Hits_neg_pres'] = hits['Hits_neg_pres'] - hits['DAP_Hits_neg_pres']
  hits = hits.fillna(0)

  hits['temp'] = hits.apply(lambda r: stats.fisher_exact([[r.DAP_Hits_pos_pres, r.noDAP_Hits_pos_pres], [r.DAP_Hits_neg_pres, r.noDAP_Hits_neg_pres]]), axis=1)
  hits[['DAP_Odds', 'DAP_pval']] = hits['temp'].apply(pd.Series)
  hits['DAP%_OnlyPos'] = hits['DAP_Hits_pos_pres']/hits['Hits_pos_pres']
  hits['DAP%'] = (hits['DAP_Hits_pos_pres'] + hits['DAP_Hits_neg_pres'])/ (hits['Hits_pos_pres'] + hits['Hits_neg_pres'])

  
  df = pd.merge(df, hits.loc[:, ['DAP_Odds', 'DAP_pval', 'DAP%', "DAP%_OnlyPos"]], how = 'left', right_index=True, left_index=True)

if FS:
  def overlap_function(x):
    for six_mer in fs_kmers:
      if six_mer in x:
        return 1
    return 0
  fs = pd.read_csv(FS, sep = "\t", header = 0, index_col=0)
  fs_kmers = np.delete(fs.columns.values, 0)
  df['FS_ovlp'] = df.index.values
  df['FS_ovlp'] = df['FS_ovlp'].apply(overlap_function)


if TFBM_CISBP:
  cisbp = pd.read_csv(TFBM_CISBP, sep='\t', header=0, index_col=0)
  min_val =  pd.DataFrame(cisbp.min(axis=1), columns = ['CISBP_TFBM_PCCD'])
  min_id = pd.DataFrame(cisbp.idxmin(axis=1), columns = ['CISBP_TFBM'])

  df = pd.merge(df, min_val, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, min_id, how = 'left', right_index=True, left_index=True)

if TFBM_DAPSeq:
  dapseq = pd.read_csv(TFBM_DAPSeq, sep='\t', header=0, index_col=0)
  min_val =  pd.DataFrame(dapseq.min(axis=1), columns = ['DAPSeq_TFBM_PCCD'])
  min_id = pd.DataFrame(dapseq.idxmin(axis=1), columns = ['DAPSeq_TFBM'])

  df = pd.merge(df, min_val, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, min_id, how = 'left', right_index=True, left_index=True)




print(df.head())
