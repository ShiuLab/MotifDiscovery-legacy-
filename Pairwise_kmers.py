"""Set up and run datasets through Random Forest

FUNCTIONS (-f):

make_pairs :=  		Make a file with all possible kmer pairs. For 6mers, should have 8,386,560 pairs. Need: -k
make_df    := 	 	Make a table with presence or absense of all kmers/kmer pairs for positive and negative genes. 
			  		Works as input for RandomForest in R. Need: -kmers, -p (fasta files), -n (fasta file), -ds, -out. 
			  			If no DNA Structure "-ds no"


PARAMETERS AVAILABLE:

-f 			:=		functions - defined above
-k 			:=		list of kmers. For make_pairs, list of singletons. For make_df, all k-mers/pairs you want in df
-p			:=		fasta files of positive examples (Will be given Class = 1)
-n 			:=		fasta files of negative examples (Will be given Class = 0)
-df 		:=		Presence or Absense dataframe (like those made with make_df).
-ds 		:=		DNA Structural information (output from /mnt/home/azodichr/scripts/...)
-pval 		:=		Default = 0.05


'''

"""

from collections import defaultdict
import sys, os
from Bio import SeqIO
import pandas as pd
import numpy as np
from scipy import stats as stats

class Kmer_pairs:

	def make_pairs(self,kmers):
		"""Make a file with all possible kmer pairs. For 6mers, should have 8,386,560 pairs"""

		km = []

		for l in open(kmers,'r'):
			k = l.strip("\n")
			km.append(k)


		pairs =[]
		for i in range(0,len(km)-1):
			for j in range(1,len(km)):
				if i+j < len(km):
					pairs.append(km[i]+ " "+km[i+j])

		print("kmers: "+str(len(km)) + ", Pairs: " + str(len(pairs)))

		out = open(kmers+"_pairs.txt",'w')
		for p in pairs:
			out.write("\n"+p)

	def make_df(self, kmers, pos, neg, ds):
		"""Make a table with presence or absense of all kmers/kmer pairs for 
		positive and negative genes. For input into randomForest"""
		n = pos.strip().split("/")[-1]
		na = n[:-7]

		km = []
		for l in open(kmers, 'r'):
			km.append(l.strip("\n"))

		##Read positive fasta files into dictionary
		genes = {}
		p = open(pos,'r')
		for seq_record in SeqIO.parse(p, 'fasta'):
			header = seq_record.id
			seq = (str(seq_record.seq))
			genes[header]=seq
		print("Positive Fasta file loaded")

		##Read neg fasta files into dictionary
		genes_neg = {}
		n = open(neg,'r')
		for seq_record in SeqIO.parse(n, 'fasta'):
			header = seq_record.id
			seq = (str(seq_record.seq))
			genes_neg[header]=seq
		print("Negative Fasta file loaded")

		dsinfo = defaultdict(list)
		if not ds.startswith('no'): 
			for l in open(ds,'r'):
				gene = l.split('\t')[0]
				DNAS = l.strip().split('\t')[1:6]
				dsinfo[gene]=DNAS


		##Mark presence or absense of motifs in each sequence in the fasta file
		m = 0
		allgenes = defaultdict(list)
		for pi in genes:
			m+=1
			templist = []
			templist.append("1")
			for ki in km:
				if ki=="DS5":
					info= '\t'.join(dsinfo[pi])
					templist.append(info)
				elif " " in ki:
					x = (len(ki)-1)/2
					out_name = na+'_'+str(x)+"paired_df.txt"
					kmer1 = ki.split(" ")[0]
					kmer2 = ki.split(" ")[1]
					seq = genes[pi]
					if kmer1 in seq and kmer2 in seq:
						templist.append("1")
					else:
						templist.append("0")
				else:
					x = len(ki)
					out_name = na+'_'+str(x)+"single_df.txt"
					seq = genes[pi]
					if ki in seq:
						templist.append("1")
					else:
						templist.append("0")
						
			allgenes[pi]=templist
			if m%25==0:
				print("Completed " + str(m) + " positive sequences")
			
		print("All Positive Examples in Dictionary")

		j = 0
		for ni in genes_neg:
			j+=1
			templist = []
			templist.append("0")
			for ki in km:
				if ki=="DS5":
					info= '\t'.join(dsinfo[ni])
					templist.append(info)
				elif " " in ki:
					kmer1 = ki.split(" ")[0]
					kmer2 = ki.split(" ")[1]
					seq = genes_neg[ni]
					if kmer1 in seq and kmer2 in seq:
						templist.append("1")
					else:
						templist.append("0")
				else:
					seq = genes_neg[ni]
					if ki in seq:
						templist.append("1")
					else:
						templist.append("0")
						
			allgenes[ni]=templist
			if j%25==0:
				print("Completed " + str(j) + " negative sequences")

		print("All Negative Examples in Dictionary")

		out = open(out_name,'w')
		if not ds.startswith('no'):
			out.write("Gene\tClass\t"+"DS1\tDS2\tDS3\tDS4\t"+'\t'.join(km))
		else:
			out.write("Gene\tClass\t"+'\t'.join(km))
		for alli in allgenes:
			out.write("\n"+alli +"\t"+ "\t".join(allgenes[alli]))
		print("Done! # genes:" + str(len(allgenes)))

	def parse(self, pos, neg, pval):
		"""Parse table based on enrichment using Fishers Exact Test, default p value is 0.05"""
		n = pos[:-4]
		out = open(n+"_FET.txt",'w')



		positives = defaultdict(list)
		for line in open(pos,'r'):
			if line.startswith("Gene"):
				header = line.strip("\n").strip("\r").split("\t")
			else:
				x = line.strip().split("\t")
				positives[x[0]]=x[2:]
		print("Positive examples loaded into dictionary")

		num_pos = len(positives)
		df_p=pd.DataFrame(positives, dtype="float",index=header[2:])
		df_p = df_p.T 				#Transpose dataframe so motif pairs are separated by column and genes by row.
		print("Positive examples loaded into panda df")
		SUMS_p = df_p.sum(0)
		print("Positive example sums calculated")

		negatives = defaultdict(list)
		for line in open(neg,'r'):
			if line.startswith("Gene"):
				header = line.strip("\n").strip("\r").split("\t")
			else:
				x = line.strip().split("\t")
				negatives[x[0]]=x[2:]
		print("Negative examples loaded into dictionary")
		num_neg = len(negatives)
		df_n=pd.DataFrame(negatives, dtype="float",index=header[2:])
		df_n = df_n.T 				#Transpose dataframe so motif pairs are separated by column and genes by row.
		print("Negative examples loaded into panda df")
		SUMS_n = df_n.sum(0)
		print("Negative example sums calculated")

		count = 0
		for i in header[2:]:
			count +=1
			oddsratio,pvalue = stats.fisher_exact([[SUMS_p[i],SUMS_n[i]],[num_pos-SUMS_p[i],num_neg-SUMS_n[i]]])
			out.write(i + "\t" + str(SUMS_p[i])+ "\t" + str(SUMS_n[i])+"\t" + str(pvalue) + "\n")
			if float(pvalue) > float(pval): ###Change back to pvalue instead of 0.5
				del df_n[i]
				del df_p[i]

			if count%5000==0:
				print("Completed " + str(count) + " kmer pairs.")

		df_p.to_csv("testout_p.csv")
		df_n.to_csv("testout_n.csv")

				
	def parse2(self, df, pval):
		"""Parse table based on enrichment using Fishers Exact Test, default p value is 0.05. Uses df with both pos and neg
		examples as the imput"""

		n = df[:-4]
		out = open(n+"_FET.txt",'w')

		num_pos = 0
		num_neg = 0
		positives = defaultdict(list)
		negatives = defaultdict(list)

		for line in open(df,'r'):
			if line.startswith("Gene"):
				header = line.strip("\n").strip("\r").split("\t")
			else:
				x = line.strip().split("\t")
				if x[1] == "1":
					positives[x[0]]=x[2:]
					num_pos+=1
				if x[1] == "0":
					negatives[x[0]]=x[2:]
					num_neg+=1
		print("Examples loaded into dictionary")

		df_p=pd.DataFrame(positives, dtype="float",index=header[2:])
		df_p = df_p.T 				#Transpose dataframe so motif pairs are separated by column and genes by row.
		print("Positive examples loaded into panda df")
		SUMS_p = df_p.sum(0)
		print("Positive example sums calculated")

		df_n=pd.DataFrame(negatives, dtype="float",index=header[2:])
		df_n = df_n.T 				#Transpose dataframe so motif pairs are separated by column and genes by row.
		print("Negative examples loaded into panda df")
		SUMS_n = df_n.sum(0)
		print("Negative example sums calculated")

		missed = []
		count = 0
		for i in header[2:]:
			count +=1
			try:
				oddsratio,pvalue = stats.fisher_exact([[SUMS_p[i],SUMS_n[i]],[num_pos-SUMS_p[i],num_neg-SUMS_n[i]]])
				out.write(i + "\t" + str(SUMS_p[i])+ "\t" + str(SUMS_n[i])+"\t" + str(pvalue) + "\n")
			except ValueError:
				missed.append(i)
			#if float(pvalue) > float(pval): ###Change back to pvalue instead of 0.5
				#del df_n[i]
				#del df_p[i]

			if count%50000==0:
				print("Completed " + str(count) + " kmer pairs.")
		print (str(len(missed))+": motif pairs were skipped because they had ValueError messages")
		#nname = n+"_neg.csv"
		#pname = n+"_pos.csv"
		#df_p.to_csv(pname)
		#df_n.to_csv(nname)



#-------------------------------------------------------------------------------
if __name__ == "__main__":
	Kmer_pairs=Kmer_pairs()
    # Print the main function help if no inputs are given.
    	pval = 0.05
	for i in range (1,len(sys.argv),2):
	    if sys.argv[i] == "-k":
	    	kmers = sys.argv[i+1]
	    elif sys.argv[i]=="-f":
	    	F = sys.argv[i+1]
	    elif sys.argv[i]=="-df":
	        df = sys.argv[i+1]
	    elif sys.argv[i]=="-p":
	        pos = sys.argv[i+1]
	    elif sys.argv[i]=="-n":
	        neg = sys.argv[i+1]
	    elif sys.argv[i]=="-ds":
	        ds = sys.argv[i+1]
	    elif sys.argv[i]=="-pval":		#default set at 0.05
	        pval = sys.argv[i+1]
	    else:
	    	print "Unknown flag: ",sys.argv[i]

	if len(sys.argv) <= 1:
		print __doc__
		sys.exit()

	if F == "make_pairs":
		if "" in [kmers]:
			print "Need files with all k-mers"
                Kmer_pairs.make_pairs(kmers)

        elif F == "make_df":
                if "" in [kmers, pos, neg]:
                        print "Need kmer list (single or pairs), pos and neg fasta files, DNA strucuture file if desired, and output directory."
                Kmer_pairs.make_df(kmers, pos, neg, ds)
    
        elif F == "parse":
                if "" in [pos, neg]:
                        print "Need presenese/absense data frame for positive and negative examples and output directory"
                Kmer_pairs.parse(pos, neg, pval)
	elif F == "parse2":
                if "" in [df]:
                        print "Need presenese/absense data frame with both positive and negative examples (& Class Column) and output directory"
                Kmer_pairs.parse2(df, pval)




