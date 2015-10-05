# MotifDiscovery
Working with RandomForest in python and R to improve motif discovery

Pairwise_kmers.py: Contains functions to make lists of paired kmers, make data frames of what genes contain those motifs, and run Fisher's Exact test to determine enrichment of those kmers/kmer pairs in the positive genes. 
RandomForest.R: Runs Random Forest on input dataframe. 10 replicates and 10 fold cross validation. 


# What you need:
•	File with all your positive examples (naming scheme will be based off the name of the positive example file, so make sure that makes sense and isn't too long)
•	File with all your negative examples
•	Fasta file with all gene promoter regions. 

If you already have fasta files of positive and neg examples, skip steps 2-3.

In HPC load:  Python3, Biopython, and SciPy

Scripts: /mnt/home/azodichr/GitHub/MotifDiscovery/

# GETTING YOUR FILES SET UP:
1. Inside directory for Pairwise experiment make directories for FASTA files and Motif Lists:
$ mkdir FastaFiles
$ mkdir MotifLists

2. Put cluster file in FastaFiles dir and get promoter sequence:
$ cd FastaFiles/
$ cp [pos_examples] .
$ python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta [promoter_sequences*] -name [pos_examples]
*For arabidopsis you can use: /mnt/home/azodichr/01_DualStress_At/TAIR10_upstream1000_Alex_20101104.mod.fa

3. Put negative example file in FastaFiles dir and get promoter sequence. 
$ cp [neg_examples] .     #For Random Forest you want a 1:1 ratio of positive and negative examples, if 73 genes in cluster, randomly select 73 negative examples.
$ python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta [promoter_sequences*] -name [neg_examples]
*For arabidopsis you can use: /mnt/home/azodichr/01_DualStress_At/TAIR10_upstream1000_Alex_20101104.mod.fa

4. Make singleton and paired kmer list. Reverse complement sequences are separated by “.”, pairs separated by a space (k = length of kmer you want):
*Using these codes you will get singleton and pair lists for 5-mers and 6-mers.
$ cd ../MotifLists
$ python Pairwise_kmers.py -f make_pairs2 –k 5
$ python Pairwise_kmers.py -f make_pairs2 –k 6


# MAKE YOUR PRESENCE/ABSENCE DATAFRAMES AND PARSE FOR ENRICHMENT
The work flow here is 1) Make df with all kmers/pairs. 2) Make list of enriched kmers/pairs. 3) Remake df with just those enriched kmers/pairs.

1. Make data frame with presence or absence of all kmer/kmer pair:
$ python Pairwise_kmers.py -f make_df –k [ListOfKmers] -p [positive fasta files] -n [negative fasta files]
*If you want to add DNA structure information to your prediction ask me. It didn't add much to my prediction...

2. If you want all the motifs move on to step 4, otherwise parse your motifs using Fisher's Exact Test:
$ python Pairwise_kmers.py -f parse_df –df [output df from step 5]
  OPTIONAL: -pval <Default is 0.05>

3. Re-make data frame with only enriched motifs:
$ python Pairwise_kmers.py -f make_df –k [output from step 6, ending in: “_sig_0.05.txt”] -p [positive fasta files] -n [negative fasta files]

4. Run Random Forest in R. If randomForest is not in your library yet, see the end of this document.
$ export R_LIBS_USER=~/R/library
$ R --vanilla --slave --args [df*] < RandomForest.R
*Can use output df from step 1 or 3

This will output two files:
.imp.txt: Open in exel, sort by "Mean Decrease Accuracy" - Make sure you shift the column headers over by one- they skip the motif name heading...
.Results.txt: Output with F-measure, stdev, sterror, and 95% confidence intervals.



Getting RandomForest onto HPC:
$ Rscript -e "install.packages(‘LIBRARY_NAME',lib='~/R/library',contriburl=contrib.url('http://cran.r-project.org/'))”
$ export R_LIBS_USER=~/R/library      #you will need to run this line every time you run RandomForest.R
$ R
$ Sys.getenv("R_LIBS_USER")
$ library(“LIBRARY_NAME")
$ q()
