"""
PURPOSE:
Given a list of kmers, a folder with mapping files, and a set of genes, list their mapping locations

Before running, set path to Miniconda:  export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH


INPUT:
  -pos_file FASTA file with positive examples


OUTPUT:
  -S



python ~/GitHub/MotifDiscovery/MapEnrichedKmers.py -k ../01_RF_Results/16_2_5_UpdatedRF/L_NNU_NNN_imp.txt -tamo /mnt/home/azodichr/01_DualStress_At/14_LogicClusters/06_MapKmers/02_mapped/ -genes ../01_RF_Results/16_2_5_UpdatedRF/00_SeqFiles/L_NNU.txtls
  """

from collections import defaultdict
import pandas as pd
import sys, os

class DAP_Seq:

  def map_specific_kmers(self, DF, PATH, POS, NEG):
    
    #Load feature matrix and save feature names 
    if isinstance(DF, str):                               #If loading df as a file
      df = pd.read_csv(DF, sep='\t', index_col = 0)   
    else:                                                 #If df is already in pandas format from previous script
      df = DF

    wanted_kmers = list(df.columns.values)
    
    pos_kmers = defaultdict(list)
    neg_kmers = defaultdict(list)

    #Make dictionary of positive example genes
    #key = kmer, value = list of all genes with that kmer & its location (GENE:Start:Stop:Strand)
    for i in os.listdir(PATH):
      if POS in i and "0.9max.out.pvalue" in i:
        k_list = []
        for l in open(os.path.join(PATH,i),'r'):
          if l[0].isdigit():
            kmer = l.strip().split("\t")[1]
          else:
            loci = ":".join(l.strip().split("\t")[0:4])
            pos_kmers[kmer].append(loci)

    #Repeat to make dictionary with negative example genes.
    for i in os.listdir(PATH):
      if NEG in i and "0.9max.out.pvalue" in i:
        k_list = []
        for l in open(os.path.join(PATH,i),'r'):
          if l[0].isdigit():
            kmer = l.strip().split("\t")[1]
          else:
            loci = ":".join(l.strip().split("\t")[0:4])
            neg_kmers[kmer].append(loci)
    print(pos_kmers)

    

    #Load bed file

    #For each kmer record the highest, lowest, and average availability amongst the genes in the pos vs. neg list. Calculate mean dif. 

    #Search bed file for positive and negative genes with kmers present

    #Output
    

  def all_DAP_sites(self, DAP, POS, NEG):
    if isinstance(DAP, str):                               #If loading df as a file
      dap = pd.read_csv(DAP, sep='\t')   
    else:                                                 #If df is already in pandas format from previous script
      dap = DAP
  
    with open(POS) as p:
      pos = p.read().splitlines()

    with open(NEG) as n:
      neg = n.read().splitlines()

    df_p = dap[dap['Gene'].isin(pos)]
    df_p.insert(1, 'Class', 1)
    num_pos = df_p.shape[0]
    df_n = dap[dap['Gene'].isin(neg)]
    df_n.insert(1, 'Class', 0)
    num_neg = df_n.shape[0]

    frames = [df_p, df_n]
    df = pd.concat(frames)

    for column in df:
      if column == "Class" or column == "Gene":
        pass
      else:
        pass
        #print(df[column])
    
    df.to_csv("all_DAP_sites_co2k4_NND_df.txt", sep="\t", index=False)



#-------------------------------------------------------------------------------

if __name__ == '__main__':
  DAP_Seq=DAP_Seq() # Print main function help if no imputs given

  for i in range (1,len(sys.argv),2):

        if sys.argv[i] == '-f':      #Directory with mapped tamo files
          F = sys.argv[i+1]
        if sys.argv[i] == '-maps':      #Directory with mapped tamo files
          PATH = sys.argv[i+1]
        if sys.argv[i] == '-pos':       #String that will identify positive genes in the maps directory 
          POS = sys.argv[i+1]
        if sys.argv[i] == '-neg':       #String that will identify negative genes in the maps directory 
          NEG = sys.argv[i+1]
        if sys.argv[i] == '-df':        #Dataframe of pos-neg examples and kmer features with presence and absense of kmers.
          DF = sys.argv[i+1]
        if sys.argv[i] == '-DAP_matrix':   #DF with all DAP TFs as columns and all genes as rows, 1 = hit in promoter
          DAP = sys.argv[i+1]

  if len(sys.argv) <= 1:
    print(__doc__)
    sys.exit()
  
  if F == "map_specific_kmers":
    if "" in [DF, PATH, POS, NEG]:
      print("Need df, path, pos, and neg")
    DAP_Seq.map_specific_kmers(DF, PATH, POS, NEG)

  elif F == "all_DAP_sites":
    if "" in [DAP, POS, NEG]:
      print("Need df, path, pos, and neg")
    DAP_Seq.all_DAP_sites(DAP, POS, NEG)



