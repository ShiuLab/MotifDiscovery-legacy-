# MotifDiscovery
Working with RandomForest in python and R to improve motif discovery

Pairwise_kmers.py: Contains functions to make lists of paired kmers, make data frames of what genes contain those motifs, and run Fisher's Exact test to determine enrichment of those kmers/kmer pairs in the positive genes. 



What you need:
•	File with all your positive examples (naming scheme will be based off the name of the positive example file, so make sure that makes sense and isn't too long)
•	File with all your negative examples
•	Fasta file with all gene promoter regions. 

If you already have fasta files of positive and neg examples, skip steps 2-3.

In HPC load:  Python3, Biopython, and SciPy

Scripts: /mnt/home/azodichr/GitHub/MotifDiscovery/

1. Inside directory for Pairwise experiment make directory for FASTA files:
$ mkdir FastaFiles

2. Put cluster file in FastaFiles dir and get promoter sequence:
$ cp ~/01_DualStress_At/00_Clusters/00_RasCat_Kmer/73.co1k3.txt FastaFiles/
$ python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta ~/01_DualStress_At/TAIR10_upstream1000_Alex_20101104.mod.fa -name FastaFiles/73.co1k3.txt

3. Select negative example files in order to have balanced experiment:
$ head -73 ~/01_DualStress_At/negativecontrol_good_1000.txt > 73.neg.txt
$ python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta ~/01_DualStress_At/TAIR10_upstream1000_Alex_20101104.mod.fa -name 73.neg.txt

4. Make singleton and paired kmer list. Reverse complement sequences are separated by “.”, pairs separated by a space (k = length of kmer you want):
$ python Pairwise_kmers.py -f make_pairs2 –k 6

5. Make data frame with presence or absence of all kmer/kmer pair:
$ python Pairwise_kmers.py -f make_df2 –k <ListOfKmers> -p <positive fasta files> -n <negative fasta files> -ds <DNA Structure info/no> 

6. If you want all the motifs move on to step 8, otherwise parse your motifs using Fisher's Exact Test:
$ python Pairwise_kmers.py -f parse_df –df <df from step 5> -pval <Default is 0.05>

7. Re-make data frame with only enriched motifs:
$ python Pairwise_kmers.py -f make_df2 –k <output from step 6, ending in: “_sig_0.05.txt” > -p <positive fasta files> -n <negative fasta files> -ds <DNA Structure info/no> 

8. Load data frame (all or fisher’s enriched) into R:
