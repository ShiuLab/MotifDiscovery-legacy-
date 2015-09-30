"""Set up and run datasets through Random Forest

FUNCTIONS (-f):

make_pairs 	:=  	Make a file with all possible kmer pairs. For 6mers, should have 8,386,560 pairs.
						Need: -k (# for k)

make_pairs2 :=  	Make pairs and singleton file taking into consideration reverse complements. 
						Need: -k (# for k)

make_df    	:= 	 	Make a table with presence or absense of all kmers/kmer pairs for positive and negative genes. 
			  		Works as input for RandomForest in R. If including DNA Structure, make sure "DS5" is in your kmer list!!!!
			  			Need: -k, -p (fasta files), -n (fasta file), -ds. 
			  				*If no DNA Structure "-ds no"
parse2		:=		Parse table based on enrichment using Fishers Exact Test, default p value is 0.05. Uses df with both pos and neg
					examples as the imput. 
						Need: -df
						Optional: -pval (Default = 0.05)
parse		:=		Old implementation of parse2, based on when there were separate df for the pos and neg examples.
						Need: -p, -n
						Optional: -pval (Default = 0.05)

PARAMETERS AVAILABLE:

-f 			:=		functions - defined above
-k 			:=		Varies by function. 
						-make_pairs: # of kmer you want
						-make_df: txt file with list of all k-mers/pairs you want in the df
-p			:=		fasta files of positive examples (Will be given Class = 1)
-n 			:=		fasta files of negative examples (Will be given Class = 0)
-df 		:=		Presence or Absense dataframe (like those made with make_df).
-ds 		:=		DNA Structural information (output from /mnt/home/azodichr/scripts/...)
-pval 		:=		Default = 0.05



"""

from collections import defaultdict
import sys, os
from Bio import SeqIO
import itertools
from Bio.Seq import Seq
#import pandas as pd
#import numpy as np
#from scipy import stats as stats

class Kmer_pairs:

	def make_pairs(self,kmers):
		"""Make a file with all possible kmer pairs. For 6mers, should have 8,386,560 pairs"""

		#Makes list of all possible kmers
		bases = ['A','T','G','C']
		km = [''.join(p) for p in itertools.product(bases, repeat=int(kmers))]

		#Make all pairwise combinations of the pairs - order does not matter and a kmer does not pair with itself
		pairs =[]
		for i in range(0,len(km)-1):
			for j in range(1,len(km)):
				if i+j < len(km):
					pairs.append(km[i]+ " "+km[i+j])

		print("kmers: "+str(len(km)) + ", Pairs: " + str(len(pairs)))

		out = open(kmers+"_pairs.txt",'w')
		for p in pairs:
			out.write("\n"+p)



	def make_pairs2(self,kmers):
		"""Make a file with all possible kmer pairs accounting for reverse complements. RC separated by '.', pairs by ' '"""

		#Makes list of all possible kmers
		bases = ['A','T','G','C']
		km = [''.join(p) for p in itertools.product(bases, repeat=int(kmers))]
		print("Possible kmers: " + str(len(km)))
		
		#Removes reverse complements so only one version present in list
		for i in km:
			s = Seq(i)
			if s.reverse_complement() in km:
				km.remove(s.reverse_complement())
		print("Kmers (reverse complements removed): " + str(len(km)))

		#Make list of all kmers with their reverse complement 
		rc_list = []
		for j in km:
			revcomp = Seq(j).reverse_complement()
			string = str(j) + "." + str(revcomp)
			rc_list.append(string)

		out = open(kmers+"mers_withRC.txt",'w')
		for k in rc_list:
			out.write("%s\n" % k)

		#Make all pairwise combinations of the pairs - order does not matter and a kmer does not pair with itself
		pairs =[]
		for i in range(0,len(rc_list)-1):
			for j in range(1,len(rc_list)):
				if i+j < len(rc_list):
					pairs.append(rc_list[i]+ " "+rc_list[i+j])

		print("Number of pairs generated: kmers: "+ str(len(pairs)))

		out2 = open(kmers+"mer_pairs_withRC.txt",'w')
		for p in pairs:
			out2.write("%s\n" % p)

	def make_df(self, kmers, pos, neg, ds):
		"""Make a table with presence or absense of all kmers/kmer pairs for positive and negative genes. 
		For input into randomForest. If inlcuding DNA Structure, include "DS5" in your kmer list"""

		#Get name for saving df, based on positive fasta file name. 
		n = pos.strip().split("/")[-1]
		na = n[:-7]

		#Put all kmers/kmer pairs into list
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

		##Read DNA Structure information into dictionary --- if -ds != no
		dsinfo = defaultdict(list)
		if not ds.startswith('no'): 
			for l in open(ds,'r'):
				gene = l.split('\t')[0]
				DNAS = l.strip().split('\t')[1:6]	#Each row has gene name and then the 5 DS principle components
				dsinfo[gene]=DNAS


		##Mark presence or absense of motifs in each sequence in the pos fasta file - for each gene go through all motifs/DNA Structure
		m = 0
		allgenes = defaultdict(list)
		for pi in genes:
			m+=1
			templist = []
			templist.append("1")
			for ki in km:
				if ki=="DS5":			#If you're including DS in dataframe, you should includ "DS5" in your motif list
					info= '\t'.join(dsinfo[pi])
					templist.append(info)	
				elif " " in ki:			#Checks to see if motif is a pair - pairs are separated by a space
					x = (len(ki)-1)/2
					out_name = na+'_'+str(x)+"paired_df.txt"
					kmer1 = ki.split(" ")[0]
					kmer2 = ki.split(" ")[1]
					seq = genes[pi]
					if kmer1 in seq and kmer2 in seq:
						templist.append("1")
					else:
						templist.append("0")
				else:					#If not DS5, and no separation by a space, assumes you're looking at singletons.
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

		##Mark presence or absense of motifs in each sequence in the neg fasta file - for each gene go through all motifs/DNA Structure

		j = 0
		for ni in genes_neg:
			j+=1
			templist = []
			templist.append("0")
			for ki in km:
				if ki=="DS5":			#If you're including DS in dataframe, you should includ "DS5" in your motif list
					info= '\t'.join(dsinfo[ni])
					templist.append(info)
				elif " " in ki:			#Checks to see if motif is a pair - pairs are separated by a space
					kmer1 = ki.split(" ")[0]
					kmer2 = ki.split(" ")[1]
					seq = genes_neg[ni]
					if kmer1 in seq and kmer2 in seq:
						templist.append("1")
					else:
						templist.append("0")
				else:					#If not DS5, and no separation by a space, assumes you're looking at singletons.
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


	def make_df2(self, kmers, pos, neg, ds):
		"""Make a table with presence or absense of all kmers/kmer pairs for positive and negative genes. 
		This version works for kmer list that accound for reverse complements
		For input into randomForest. If inlcuding DNA Structure, include "DS5" in your kmer list"""

		#Get name for saving df, based on positive fasta file name. 
		n = pos.strip().split("/")[-1]
		na = n[:-7]

		#Put all kmers/kmer pairs into list
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

		##Read DNA Structure information into dictionary --- if -ds != no
		dsinfo = defaultdict(list)
		if not ds.startswith('no'): 
			for l in open(ds,'r'):
				gene = l.split('\t')[0]
				DNAS = l.strip().split('\t')[1:6]	#Each row has gene name and then the 5 DS principle components
				dsinfo[gene]=DNAS


		##Mark presence or absense of motifs in each sequence in the pos fasta file - for each gene go through all motifs/DNA Structure
		m = 0
		allgenes = defaultdict(list)
		for pi in genes:
			m+=1
			templist = []
			templist.append("1")
			for ki in km:
				if ki=="DS5":			#If you're including DS in dataframe, you should includ "DS5" in your motif list
					info= '\t'.join(dsinfo[pi])
					templist.append(info)	
				elif " " in ki:			#Checks to see if motif is a pair - pairs are separated by a space
					x = (len(ki)-1)/2
					out_name = na+'_'+str(x)+"paired_df.txt"
					kmer1 = ki.split(" ")[0]
					kmer2 = ki.split(" ")[1]
					seq = genes[pi]
					if kmer1 in seq and kmer2 in seq:
						templist.append("1")
					else:
						templist.append("0")
				else:					#If not DS5, and no separation by a space, assumes you're looking at singletons.
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

		##Mark presence or absense of motifs in each sequence in the neg fasta file - for each gene go through all motifs/DNA Structure

		j = 0
		for ni in genes_neg:
			j+=1
			templist = []
			templist.append("0")
			for ki in km:
				if ki=="DS5":			#If you're including DS in dataframe, you should includ "DS5" in your motif list
					info= '\t'.join(dsinfo[ni])
					templist.append(info)
				elif " " in ki:			#Checks to see if motif is a pair - pairs are separated by a space
					kmer1 = ki.split(" ")[0]
					kmer2 = ki.split(" ")[1]
					seq = genes_neg[ni]
					if kmer1 in seq and kmer2 in seq:
						templist.append("1")
					else:
						templist.append("0")
				else:					#If not DS5, and no separation by a space, assumes you're looking at singletons.
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
        if F == "make_pairs2":
		if "" in [kmers]:
			print "Need files with all k-mers"
                Kmer_pairs.make_pairs2(kmers)

        elif F == "make_df":
                if "" in [kmers, pos, neg]:
                        print "Need kmer list (single or pairs), pos and neg fasta files, DNA strucuture file if desired, and output directory."
                Kmer_pairs.make_df(kmers, pos, neg, ds)
        
        elif F == "make_df2":
                if "" in [kmers, pos, neg]:
                        print "Need kmer list (single or pairs), pos and neg fasta files, DNA strucuture file if desired, and output directory."
                Kmer_pairs.make_df2(kmers, pos, neg, ds)
    
        elif F == "parse":
                if "" in [pos, neg]:
                        print "Need presenese/absense data frame for positive and negative examples and output directory"
                Kmer_pairs.parse(pos, neg, pval)
	elif F == "parse2":
                if "" in [df]:
                        print "Need presenese/absense data frame with both positive and negative examples (& Class Column) and output directory"
                Kmer_pairs.parse2(df, pval)




