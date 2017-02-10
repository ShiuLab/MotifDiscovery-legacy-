"""

Merge kmer information into one database. 
Inputs:

python ~/GitHub/MotifDiscovery/pCRE_criteria.py -df ../NDD_Ras_FET01_df_p0.01.txt -WSC ../02_kmer_mapping/01_map_WithinSpeciesDiv/all_kmers.bed.SeqStats -enrich NDD_Ras_FET01_FETresults.txt -imp_RF NDD_Ras_FET01_RF_imp -imp_SVC NDD_Ras_FET01_SVC_imp -map ../02_kmer_mapping/all_kmers.gff -DHS ../02_kmer_mapping/all_kmers.gff.DHS_all_ovrlp -CNS ../02_kmer_mapping/all_kmers.gff.CNS_between_ovrlp  -FS NDD_Ras_df.txt_decisiontree_50 -TFBM_DAPSeq ../03_KnownTFBM_Similarity/02_DAP_Seq/kmers.txt.tamo-DAP_motifs.txt.tm_mod.txt -TFBM_CISBP ../03_KnownTFBM_Similarity/01_CIS_BP/kmers.txt.tamo-Athaliana_TFBM_v1.01.tm.index.direct.index.tm_mod.txt -Histone ../02_kmer_mapping/02_Histone/histone_overlap_files.txt

##### Input #####
Required:
-df                 Dataframe with Col 1 = genes, Col 2 = Class, Col 3-... = pCREs

Optional:
-enrich             FET results from pCRE finding pipeline (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/01_Kmer_Features/NDD_Ras_FET01_FETresults.txt)
-imp_RF             RF Machine learning importance file 
-imp_SVC            SVC Machine learning importance file
-map                GFF file of all pCREs of interest (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/02_kmer_mapping/all_kmers.gff)
-DHS                Also required (-map)
-DAP                Also required (-map)
-CNS                Also required (-map)
-WSC                Also required (-map)
-Histone            Txt file with histone overlap info (Col 1 = histone marker name, Col 2 = activator/repressor, Col 3 = full path to overlap file)
                    See: /mnt/home/azodichr/01_CombinedStress/Prasch_HD/07_6clusters/02_kmer_mapping/02_Histone/histone_overlap_files.txt
-TFBM_CISBP         PCCD between kmers and CIS-BP motifs (mjliu/kmer_5/Athaliana_TFBM_v1.01.tm.index.direct.index.tm) using output from Ming's pipeline (~mjliu/kmer_5/kmer_5.sh)
-TFBM_DAPSeq        PCCD between kmers and DAP-Seq motifs (ShiuLab/14_DAPseq/PWM_to_tamo/DAP_motifs.txt.tm) using output from Ming's pipeline (~mjliu/kmer_5/kmer_5.sh)
-FS                 pCREs selected by decision tree feature selection rather than enrichement.


##### Possible Output Columns #####
pvalue              Enrichement Score 
percent_positives   Support score - i.e. Percent of cluster genes that contain the kmer in the promoter
Imp_RF              Importance score in RF predictive models
Imp_SVC             Importance score in SVC predictive models
Mean Copy#_pos_pres Average copy number in promoter region of cluster genes where kmer is present
Mean_Copy#_neg_pres Average copy number in promoter region of non-cluster genes where kmer is present
Copy#_Tstat         T-statistic from Welch's T-test. + = greater in pos, - = greater in negative 
Copy#_pval          Significance (2-tailed)
Mean_WSC_pos        Average within species nucleotide diverity of kmer in positive example genes
Mean_WSC_neg        Average within species nucleotide diverity of kmer in negative example genes
WSC_tstat           T-statistic from Welch's T-test. + = greater in pos, - = greater in negative 
WSC_pval            Significance (2-tailed)
DAP_odds            Fisher's Exact Test comparing pCREs overlapping with DAP-Deq peaks (200bp each) in cluster and non-cluster genes
DAP_pval            Significance score for DAP FET
DHS_Odds            Fisher's Exact Test comparing pCREs overlapping with DHS sites in cluster and non-cluster genes
DHS_pval            Significance score for DHS FET
CNS_Odds            Fisher's Exact Test comparing pCREs overlapping with CNS sites in cluster and non-cluster genes
CNS_pval            Significance score for CNS FET
FS_ovlp             Does the k-mer overlap with one of the 6-mers identified to predict the cluster using feature selection
CISBP_TFBM          Best match known TFBM from the CIS-BP Database
CISBP_TFBM_PCCD     Similarity (PCC-Distance - lower = better match) to best match known TFBM from the CIS-BP Database
DAPSeq_TFBM         Best match known TFBM from DAP-Seq Database
DAPSeq_TFBM_PCCD    Similarity (PCC-Distance - lower = better match) to best match known TFBM from DAP-Seq Database


"""
import sys
import numpy as np
np.set_printoptions(threshold=np.inf)
import pandas as pd
from scipy import stats

DF = ENRICH = IMP_RF = IMP_SVC = MAP = DHS = CNS = WSC = TFBM = DAP = FS = TFBM_CISBP = TFBM_DAPSeq = HIST = False


def ovrp_enriched(e_df, hits):
  """ Function to calculate enrichement of overlap with feature in cluster compared to non-cluster genes 

  Fisher's Exact Test
                                 Overlap with Feature   |  No overlap with feature
    pCRE in Cluster gene      |       t1                            t2
    pCRE in Non-Cluster gene  |       t3                            t4

"""
  subset = e_df[e_df['kmer'].isin(kmers)]
  subset = subset[subset.details.str.contains("overlap=NA") == False]
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]

  count_pos = subset_pos.groupby(['kmer','gene']).size()
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['t1'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['t3'])

  hits = pd.merge(hits, grouped_pos, how = 'left', right_index=True, left_index=True)
  hits = pd.merge(hits, grouped_neg, how = 'left', right_index=True, left_index=True)

  hits['t2'] = hits['Hits_pos_pres'] - hits['t1']
  hits['t4'] = hits['Hits_neg_pres'] - hits['t3']
  hits = hits.fillna(0)

  hits['temp'] = hits.apply(lambda r: stats.fisher_exact([[r.t1, r.t2], [r.t3, r.t4]]), axis=1)
  hits[['Odds', 'pval']] = hits['temp'].apply(pd.Series)

  return hits



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
      if sys.argv[i] == "-Histone":
        HIST = sys.argv[i+1]
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
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).mean(), columns = ['AvCopy#_pos_pres'])
  num_hits_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['Hits_pos_pres'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).mean(), columns = ['AvCopy#_neg_pres'])
  num_hits_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['Hits_neg_pres'])

  df = pd.merge(df, grouped_pos, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, grouped_neg, how = 'left', right_index=True, left_index=True)

  # Total number of kmer hits in promoters of pos and neg genes - usefor for fisher's enrichement tests.
  hits = pd.merge(num_hits_pos, num_hits_neg, how = 'left', right_index=True, left_index=True)

  # Perform Welch's T-test to look for statistical difference in pCRE copy number between pos and neg genes (non-normal dist)
  df_temp = pd.DataFrame(columns = ['Copy#_Tstat', 'Copy#_pval'], index = kmers)
  for k in kmers:
    val1 = count_pos.loc[k]
    val2 = count_neg.loc[k]
    tstat, pval = stats.ttest_ind(val1, val2, nan_policy = 'omit', equal_var = False)
    df_temp.loc[k] = [tstat, pval]

  df = pd.merge(df, df_temp, how = 'left', right_index=True, left_index=True )

if DAP:
  feat = pd.read_csv(DAP, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  feat_hits = ovrp_enriched(feat, hits)
  feat_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', "DAP_odds", "DAP_pval"]
  df = pd.merge(df, feat_hits.loc[:, ['DAP_odds', "DAP_pval"]], how = 'left', right_index=True, left_index=True)

if DHS:
  feat = pd.read_csv(DHS, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  feat_hits = ovrp_enriched(feat, hits)
  feat_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', "DHS_odds", "DHS_pval"]
  df = pd.merge(df, feat_hits.loc[:, ['DHS_odds', "DHS_pval"]], how = 'left', right_index=True, left_index=True)

if CNS:
  feat = pd.read_csv(CNS, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  feat_hits = ovrp_enriched(feat, hits)
  feat_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', "CNS_odds", "CNS_pval"]
  df = pd.merge(df, feat_hits.loc[:, ['CNS_odds', "CNS_pval"]], how = 'left', right_index=True, left_index=True)

if WSC:
  wsc = pd.read_csv(WSC, sep='\t', header = 0)
  wsc['gene_and_region'], wsc['region'], wsc['kmer'] = wsc['Region'].str.split('|').str
  wsc['chr'], wsc['start'], wsc['stop'], wsc['gene'] = wsc['gene_and_region'].str.split('_').str
  wsc['NtDiversity'] = wsc['NtDiversity'].replace(to_replace = 'NoDiffSites', value = 0)
  subset = wsc[wsc['kmer'].isin(kmers)]

  # Get just pos/neg example genes and convert NtDiv to numeric
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_pos['NtDiversity'] = pd.to_numeric(subset_pos['NtDiversity'], errors='coerce')
  subset_neg = subset[subset['gene'].isin(neg_genes)]
  subset_neg['NtDiversity'] = pd.to_numeric(subset_neg['NtDiversity'], errors='coerce')

  # Calculate means NtDiv for each kmer for genes in pos and neg datasets
  av_pos_mean = pd.DataFrame(subset_pos.groupby(['kmer'])['NtDiversity'].mean())
  av_pos_mean.columns = ['Mean_WSC_pos']
  av_neg_mean = pd.DataFrame(subset_neg.groupby(['kmer'])['NtDiversity'].mean())
  av_neg_mean.columns = ['Mean_WSC_neg']

  df = pd.merge(df, av_pos_mean, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, av_neg_mean, how = 'left', right_index=True, left_index=True)

  # Run Welch's T-test to determine if NtDiv of a kmer is statistically different between pos and neg genes.
  df_temp = pd.DataFrame(columns = ['WSC_tstat', 'WSC_pval'], index = kmers)
  for k in kmers:
    val1 = subset_pos[subset_pos['kmer'] == k]
    val2 = subset_neg[subset_neg['kmer'] == k]
    tstat, pval = stats.ttest_ind(val1['NtDiversity'], val2['NtDiversity'], nan_policy = 'omit', equal_var = False)
    df_temp.loc[k] = [tstat, pval]

  df = pd.merge(df, df_temp, how = 'left', right_index=True, left_index=True )

if HIST:
  with open(HIST, 'r') as hist_files:
    for l in hist_files:
      hist_name, act_repress, hist_file = l.strip().split('\t')
      his = pd.read_csv(hist_file, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
      his_hits = ovrp_enriched(his, hits)
      odds = hist_name + "_" + act_repress[0] + "_Odds"
      pval = hist_name + "_" + act_repress[0] + "_pval"
      his_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', odds, pval]
      df = pd.merge(df, his_hits.loc[:, [odds, pval]], how = 'left', right_index=True, left_index=True)

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

name = DF + '_pCRE_Criteria'
df.to_csv(name, sep = '\t')