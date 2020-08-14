'''Functions that are used to merge motifs in a TAMO file.

This script contains a set of functions that are used to merge motifs based on
there PCC score. The functions are based on scripts written by Cheng Zou, and
are insterted here virtually unchanged.

#need to purge then load modules first or using qsub_slurm.py
module purge 
module load TAMO/1.0
module load Python/2.7.10 

Usage:
python pcc_merge_CC.py [Function] -w working directory -t TAMO_file
     
Functions:
    check_jobs: You can run this script after you have run all of the PCC 
    matrix creation jobs. It creates a failed_jobs_cc, which you can rerun by 
    using qsub_hpc.py. The funciton will print the errors for the failed jobs. 
    Requires: -w
    
    check_outputs: Checks that all the proper outputs from a distance matrix
    run were created. More direct method then check_jobs
    Requires: -c (the command line file)

    check_square_matrix: Once the combined matrix is created, you can use this 
    function to check the matrix was formed correctly. Requires: -d, -w
        
    check_top_files: This function checks if each cluster that was supposed to 
    be merged has a .TOP output file. The script reports which jobs the cluster
    belongs to. Requires: -w
    
    combine_distance_matrix: After all jobs have run successfully, you can
    use this function to combine all of the distance matricies that were
    created. Requires: -t, -w
    
    create_cc: Splits the TAMO file, and creates command lines for generating 
    PCC score matricies. Requires: -w, -t
    
    create_cc_2: Splits two TAMO files into sub TAMO files, and then runs the
    pcc distance calculation as a series of jobs. Requires: -w, -t, -t2

    merge_runs: Once the PCC matrix is combined and checked, this function will
    create a tree for the motifs, cut the tree at the spcified height and then 
    merge the motifs bellow this height. Requires: -t, -w. Options: -h, -dist, 
    -ancestor, -target, -genome
    
    merge_runs_ancestor and merge_runs_no_ancestor: Used to run individual
    jobs from merge_runs_cc.

    merge_runs_cc: Divides the merge runs process into a series of jobs
    that you can run using qsub_hpc.py. Requires the same input as
    merge_runs
    
    merge_runs_cleanup: After running merge_runs jobs, you can use this to 
    combine all of the TOP files and the single cluster motifs. It also removes
    all of the unnecessary files from generated in the merging process. 
    Requires: -t, -w
    
    run_UPGMA: This script submits a job that creates the UPGMA tree from the 
    PCC score matrix. You need to run this script before running the merge_runs
    script. Requires: -t, -w. Options: -h, -wall, -m
    
Parameters:
    -ancestor := consider estimated ancestors of motifs(1) or not(0), 
    default is 1.
    
    -d := Merged PCC matrix
    
    -genome := As default, it is fasta file which job file on the array and 
    with promoter seqences annotated in TAIR8
    
    -h := height to cut the tree, Be aware of when the height is equal to 0.05,
    the distance between two taxon is 0.10
    
    -mem  := memory required to run UPGMA. Default is 256 gb.
    
    -t := TAMO file that you wish to merge
    
    -t2 := Second TAMO file for create_cc_2
    
    -target := the list of target genes that you are interested, as default it 
    is all stress responsive genes q < 0.05,2 fold change.
    
    -w := Working directory for the merging process
    
    -wall := walltime for running UPGMA. Default is 124 minutes.    
'''

import os
import sys
import math
import re

from TAMO import MotifTools
from TAMO.MotifTools import Motif, save_motifs


def TAMO_split(TAMO_file, motifs_per_file = 190):

    '''This function splits a TAMO into smaller files for create_cc'''
    ml=MotifTools.txt2motifs(TAMO_file)
    total=len(ml)/int(motifs_per_file) # Total number of TAMOs to generate
    by = motifs_per_file
    for i in range(total):
        print i
        print i*by+by,TAMO_file+'_n%s' % i
        save_motifs(ml[i*by:i*by+by],TAMO_file+'_n%s' % i)
    print     total*by,len(ml),TAMO_file+'_n%s' % (total)
    save_motifs(ml[total*by:len(ml)],TAMO_file+'_n%s' % (total))
    return(total)

def create_cc_for_two(wdir, TAMO_file_1, TAMO_file_2):

    '''Sets up pcc distance calculation commands for two TAMO files

    This command allow you to compare two different files. It divides each TAMO
    and then does the PCC calculation on the sub files.
    '''
    # Move to the working directory
    TAMO_file_1 = os.path.abspath(TAMO_file_1)
    TAMO_file_2 = os.path.abspath(TAMO_file_2)
    os.system("cd %s" % wdir)
    os.chdir(wdir)
    script_dir = os.path.abspath(__file__) # path to pcc_merge_CC.py script
    script_dir = '/'.join(script_dir.split('/')[:-1]) # path to script directory
    # Split the original TAMO files and get then number of split TAMO files
    n_split_1 = TAMO_split(TAMO_file_1, motifs_per_file = 100)
    n_split_2 = TAMO_split(TAMO_file_2, motifs_per_file = 100)    
    # Open the command line output file.
    oup=open("runcc","w")    
    for i in range(n_split_1+1):
        for j in range(n_split_2+1):
            # This command creates a matrix for comparing the ith TAMO to the 
            # jth TAMO.
            oup.write("python %s/3.calculate_disctance_2_file.py -i %s_n%s -j  %s_n%s --dfunc pccrange\n" % (script_dir,TAMO_file_1,i,TAMO_file_2,j))
    oup.close()

def create_cc(wdir, TAMO_file, dfunc):
    '''Spilts the TAMO file into sub files and creates PCC calculation commands.
    
    If your TAMO file is small, with less than 380 motifs, the function creates 
    a single command that you can use use to to merge the Whole TAMO file without
    the need to combine sub matricies.
    '''
    # Move to the working directory
    os.system("cd %s" % wdir)
    os.chdir(wdir)
    script_dir = os.path.abspath(__file__) # path to pcc_merge_CC.py script
    script_dir = '/'.join(script_dir.split('/')[:-1]) # path to script directory
    # Split the original TAMO file and get then number of split TAMO files
    n_split = TAMO_split(TAMO_file, motifs_per_file = 100)        
    # Open the command line output file.
    oup=open("runcc","w")
    # If the n_split is less than 2, the function gives you one command line
    # that you can submit to calculate the PCC distance for the whole TAMO file.
    # You need to skip the combine_distance_matrix step.
    if n_split < 2:
        print '\n---------------------------------------------------------------'
        print 'It looks like your TAMO file is fairly small. runcc has only one command line'
        print 'which can be submited with qsub_hpc.py. This command calculates the pcc' 
        print 'distance matrix for the entire TAMO file. It may require a walltime' 
        print 'of 60 minutes to run the whole thing.\n'
        print 'You can skip the steps that require combine_distance_matrix and'
        print 'check_square_matrix.\n'
        print 'Procede directly to the run_UPGMA step once the distance matrix is'
        print 'created.'
        print '----------------------------------------------------------------\n'    
        oup.write("module load TAMO; module load SciPy; python %s/3.calculate_distance_matrix.py \
-i %s/%s --dfunc pccrange\n" % (script_dir, wdir, TAMO_file))
    # Otherwise, it creates the command lines for a split TAMO file.    
    else:
        for i in range(n_split+1):
            # This command compares a TAMO file to itself
            oup.write("module load TAMO; module load SciPy; python %s/3.calculate_distance_matrix.py \
            -i  %s/%s_n%s --dfunc %s\n" % (script_dir, wdir,TAMO_file,i,dfunc))
            for j in range(i+1,n_split+1):
                # This command creates a matrix for comparing the ith TAMO to the 
                # jth TAMO.
                oup.write("module load TAMO; module load SciPy; python %s/3.calculate_disctance_2_file.py -i %s/%s_n%s -j  %s/%s_n%s --dfunc %s\n" % (script_dir, wdir,TAMO_file,i,wdir,TAMO_file,j, dfunc))
        oup.close()

def check_jobs(job_dir):
    '''Checks error files of jobs and creates a cc file for failed jobs.
    
    This function checks the error files of the completed files and checks if 
    they were killed. If they were killed it adds them to the failed command 
    line file.
    '''
    job_dir_list = os.listdir(job_dir)
    output = open('failed_jobs_cc', 'w')
    for file_name in job_dir_list:
        file_name = job_dir.rstrip('/') + '/' + file_name
        # Check if the file is an error file.
        if '.sh.e' in file_name:        
            success = True
            error_file = open(file_name, 'r')            
            # This loop ckecks for errors in the error file.
            for line in error_file:                
                if 'killed:' in line:
                    print file_name, line
                    success = False
                elif 'Traceback' in line:
                    print file_name, line
                    success = False
                elif 'Error' in line:
                    print file_name, line
                    success = False                    
            if success == False:
                job_file_name = file_name.split('.e')[0]
                job_file = open(job_file_name)
                lines = job_file.readlines()
                output.write(lines[-1])    
    output.close()

def is_good_file(fpath):
    '''Checks that a file exists and is not empty.
    '''
    return os.path.exists(fpath) and os.stat(fpath)[6] > 0

def check_outputs(wdir, cc_file):
    '''Checks that all job outputs are created and complete.
    
    This script is more direct than check_jobs because it focuses solely on the
    output files.
    '''
    ith_file_search = re.compile(r'-i\s*(.*)\s')
    jth_file_search = re.compile(r'-j\s*(.*)\s')
    dir_list = os.listdir(wdir)
    output = open('incomplete_jobs_cc', 'w')
    for line in open(cc_file, 'r'):
        #TAMO_1 = ith_file_search.match(line).groups(1)
        tab = filter(lambda x: x != '', line.strip().split(' '))
        #TAMO_1 = line.strip().split('-i')[1].strip(' ').split(' ')[0]
        for i, item in enumerate(tab):
            if item == '-i':
                TAMO_1 = tab[i+1]
            elif item == '-j':
                TAMO_2 = tab[i+1]
        if not '-j' in line:
            out_file = TAMO_1 + '.dm'
        elif '-j' in line:
            #TAMO_2 = jth_file_search.match(line).groups(1)
            #TAMO_2 = line.strip().split('-j')[1].strip(' ').split(' ')[0]
            out_file = TAMO_1 + '-' + TAMO_2.split('/')[-1] + '.dm'
        
        # Check that the proper output was generated
        if not is_good_file(out_file):
            output.write(line)
            print os.stat(out_file)[6]
            print out_file 
    output.close()
    

def line_count(file):
    '''Line count function for combine_distance_matrix'''
    count=-1
    for count,line in  enumerate(open(file)):pass
    return count+1

def remove_double_tabs(matrix):
    '''Removes double tabs in a matrix
    '''
    matrix_file = open(matrix, 'r')
    output = open(matrix + '2', 'w')
    for line in matrix_file:
        output.write(line.replace('\t\t', '\t'))
    matrix_file.close()
    output.close()

def combine_distance_matrix(wdir, TAMO_file):    
    '''Combines the PCC score matricies and outputs them as a single matrix.
    
    Originaly written by Cheng Zou, and converted to a function by Alex Seddon.
    '''    
    ml      = MotifTools.txt2motifs(TAMO_file)    
    n_split = len(ml)/100    
    ##
    # Change to the working directory.
    os.system("cd %s" % wdir)
    os.chdir(wdir)
    #
    ##    
    ##
    # The following loop keeps counts the number of lines in the each of the
    # PCC matricies for a comparison of a TAMO file with itself.    
    lendic={} # Dictionary with the length of PCC matricies.    
    for i in range(n_split + 1):
        lendic[i]= line_count("%s_n%s.dm" %  (TAMO_file,i))
    print lendic
    #
    ##    
    ##
    # This loop creates files with blanks. The files are used to ensure that 
    # the PCC-distance matrix is square. The blank files will be created to take
    # the place of files that would have been left blank 
    for i in range(n_split+1):
        for j in range(0,i):        
            # open the file to add blanks
            oup=open("%s_n%s-%s_n%s.dm" % (TAMO_file,i,TAMO_file,j),"w")
            print lendic[j],lendic[i]
            list=[]        
            # Add a number of "-" to the list equal to the number of lines in
            # the self comparison files.
            for y in range(lendic[j]):
                list.append("-")
            for x in range(lendic[i]):
                oup.write("%s\n" % "\t".join(list))
            oup.close()
    #
    ##
    
    ##
    # Creates a copy of the self comparison file so that it can be easily picked
    # out by the function.
    for i in range(n_split + 1):
        os.system("cp %s_n%s.dm %s_n%s-%s_n%s.dm" %  (TAMO_file,i,TAMO_file,i,TAMO_file,i))
    #
    ##

    ##
    # This loop will look at each 
    for i in range(n_split+1):
        com = "paste "
        for j in range(n_split+1):
            com += "%s_n%s-%s_n%s.dm " % (TAMO_file,i,TAMO_file,j)
        com += "> distance_%s" % i
        print com
        os.system(com)

    com="cat "
    for i in range(n_split+1):
        com+="distance_%s " % i
    com+= "> %s.dm" % TAMO_file
    
    print com
    # Concatonate all the matricies
    os.system(com)
    # My embarisingly ad hoc way of removing double tabs
    remove_double_tabs("%s.dm" % TAMO_file)
    #os.system("mv %s.dm2 %s.dm" % (TAMO_file, TAMO_file))
    

def check_square_matrix(wdir, PCC_matrix):

    '''This function checks the combined PCC score matrix. 
    
    It was written by Cheng Zou and converted to a function by Alex Seddon.
    '''

    os.system("cd %s" % wdir)
    os.chdir(wdir)

    os.system("cp %s distance_matrix.dm" % PCC_matrix)
    oup=open("%s.dm2" % PCC_matrix ,"w")
    n=0
    for line in open("distance_matrix.dm","r"):
        L=line.strip().split("\t")
        nL=[]
        for i in range(len(L)):
            if L[i].startswith("0") or L[i].startswith("1") or L[i].startswith("-"):
                nL.append(L[i])
            elif L[i]=="nan":
                    nL.append("1")
            else:
                print n,i,L[i-1],L[i],L[i+1]
        if len(nL)!=len(L):
            print "WRONG",len(nL),len(L)
        n+=1
        oup.write("%s\n" % "\t".join(nL))
    oup.close()

def combine_distance_matrix_for_2(wdir, TAMO_file_1, TAMO_file_2):
    
    '''Combines matricies made from two TAMO files.
    
    This script is used to create the final matrix after all jobs from 
    create_cc_for_2 are complete.
    '''
    
    ml_1      = MotifTools.txt2motifs(TAMO_file_1)
    ml_2      = MotifTools.txt2motifs(TAMO_file_2)
    
    n_split_1 = len(ml_1)/100
    n_split_2 = len(ml_2)/100
    
    print n_split_1, len(ml_1)
    print n_split_2
    
    
    # Change to the working directory.
    os.system("cd %s" % wdir)
    os.chdir(wdir)

    # This loop will paste together matricies
    for i in range(n_split_1 + 1):
        com = "paste "
        for j in range(n_split_2 + 1):
            com += "%s_n%s-%s_n%s.dm " % (TAMO_file_1,i,TAMO_file_2,j)
        com += "> distance_%s" % i
        print com
        os.system(com)

    # 
    com="cat "
    for i in range(n_split_1 + 1):
        com+="distance_%s " % i
    com+= "> %s-%s.dm" % (TAMO_file_1, TAMO_file_2)

    print com
    os.system(com)

def parse_out_STAMP(wdir,TAMO_file,i):
    
    '''This function parses the output from STAMP
    
    Originally written by Cheng Zou and converted to a funciton by Alex Seddon.
    '''
    oup=open("%s/%s_sub_%s.tm.tf_SToutFBP.txt.mod" % (wdir,TAMO_file,i),"w")
    for line in open("%s/%s_sub_%s.tm.tf_SToutFBP.txt" % (wdir,TAMO_file,i),"r"):
        if line.startswith("DE"):
            oup.write("DE\t%s-%s\n" % (TAMO_file,i))
        else:
            oup.write("%s\n" % line.strip())
    oup.close()

def parse_out_STAMP_job(TAMO_file,i):
    
    '''This function parses the output from STAMP. For use with job functions.
    
    Originally written by Cheng Zou and converted to a funciton by Alex Seddon.
    '''
    oup=open("%s_sub_%s.tm.tf_SToutFBP.txt.mod" % (TAMO_file,i),"w")
    for line in open("%s_sub_%s.tm.tf_SToutFBP.txt" % (TAMO_file,i),"r"):
        if line.startswith("DE"):
            oup.write("DE\t%s-%s\n" % (TAMO_file,i))
        else:
            oup.write("%s\n" % line.strip())
    oup.close()
#
##

##
# Help 
def help_merge_runs():
    print "-t    specify the TAMO_fileut tamo file\n"
    print "-w    specify the working directory\n"
    print "-h    height to cut the tree, Be aware of when the height is equal \
        to 0.05,the distance between two taxon is 0.10"
    print "-dist    you have already calculated distance matrix (1) or not (0)\n, \
        default is 1. if you have the matrix, it has to be named as tamo file.dm"
    print "-ancestor    consider estimated ancestors of motifs(1) or not(0), default is 1 \n"
    print "-target    the list of target genes that you are interested,\
        as default it is all stress responsive genes q<0.05,2 fold\n"
    print "-genome    As default, it is fasta file containing all the gene\
        on the array and with promoter seqences annotated in TAIR8\n"

    sys.exit(0)


def tamo2tf(TAMO_file):
    '''Converts TAMO files to the TRANSFAC format
    '''
    
    ml=MotifTools.txt2motifs(TAMO_file)
    TAMO_file_name=TAMO_file.split("/")[-1]
    ACGT=["A","C","G","T"]
    n=1
    oup=open("%s.tf" % (TAMO_file),"w")
    for m in ml:
        if m.source=="":
            oup.write("DE\t%s_%s\t%s_%s\n" % (TAMO_file_name,n,TAMO_file_name,n))
        else:
            oup.write("DE\t%s\t%s\n" % (m.source,m.source))
        count=0
        #print m.source
        for i in range(m.width):
            oup.write("%s\t" % count)
            for letter in ACGT:
                if m.logP:
                    Pij = pow(2.0, m.logP[i][letter])
                    oup.write("%s\t" % int(Pij*100))
            oup.write("\n")
            count+=1
        oup.write("XX\n")
        n+=1    
    oup.close()

##
# These two functions are used to convert TRANSFAC files to TAMO files.
def parse_block(name,block):
    mat=[]
    ACGT={"A":1,"C":2,"G":3,"T":4}
    for i in block:
        L=i.strip().split()
        D = {'A': 0, 'C': 0, 'T':0, 'G':0}
        for j in ACGT.keys():
            D[j]=float(L[ACGT[j]])
        mat.append(D)
    m=MotifTools.Motif_from_counts(mat)
    m.source=name
    #print m._print_p()
    return m

def tf2tamo(tf_file):

    inp=open(tf_file,"r")
    line =inp.readline()

    motifs=[]
    while 1:
        if not line:
            break
        if line.startswith("DE"):
            name=line.strip().split()[1]
            block=[]
            line =inp.readline()
            while not line.startswith("XX"):
                if not line :
                    break
                block.append(line)
                line =inp.readline()
            motifs.append(parse_block(name,block))
        else:
            line =inp.readline()

    save_motifs(motifs, '%s.tm' % (tf_file))
#
##


def parse_out_pcs(TAMO_file):

    '''Parses the output from the MotifMetrics.py'''
    
    from operator import itemgetter

    ##    
    # return the motif name with smallest p-value. if it is larger than 0.001,
    # return None.
    dic={}
    for line in open(TAMO_file):
        # Skips the header
        if not line.startswith("#"):
        
            # Splits the output line
            L=line.strip().split()
            
            # Gets the p-value of the motifs. 
            for i in range(len(L)):
                if L[i]=="p:":
                    dic[L[-3]]=-math.log(float(L[i+1]),10)
    
    # 
    z = sorted(dic.items(), key=itemgetter(1))
    #print "outz",z
    if z[-1][1] <2:
        return "None"
    else:
        return z[-1][0] 
#
##

def trim_motif(TAMO_file, cut = 0.4):

    '''Trims the motifs in TAMO_file, eliminating low-information flanks.'''
    
    testmotifs = MotifTools.load(TAMO_file)
    file=TAMO_file+"_"+str(cut)+".trim"

    new_mlist=[]
    for motif in testmotifs:
        m = motif.trimmed(cut)
        new_mlist.append(m)
    save_motifs(new_mlist,file)


def pick_chunk_score(wdir, TAMO_file, target, genome):

    '''Trims and returns the top motif in a cluster.
    
    This script takes in the TAMO file from the motifs in a single cluster. It
    trims the low-information ends from each motifs. It then indentifies the
    motif that is most significantly represented in the target genes in your
    genome. If no motif is significantly represented, then a blank top motif
    file is created.
    '''
    os.system("cd %s" % wdir)
    os.chdir(wdir)

    script_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1]) # path to pcc_merge_CC.py script
    
    ##
    # step 1 trim tamo to eliminate low information flanking sequence
    trim_motif(TAMO_file, 0.1)
    
    ##
    # step 2 Group Specificity Score" from the Church lab
    # python MotifMetrics.py [Genes of interest] -genome [FASTA of promoter sequence] -t [Trimmed TAMO of cluster motifs] 
    # MotifMetrics.py checks if the motifs appear disproportionatly to the 
    # targets compared to the rest of the genes.
    os.system("python %s/MotifMetrics.py %s -genome %s -t %s_0.1.trim -spec > %s_0.1.trim_Cout" % (script_dir,target,genome,TAMO_file,TAMO_file))

    ##
    # Gets the motif that is most significantly represented in your target genes
    # Returns "None" if none of the motifs has a p-value above 0.001.
    topm=parse_out_pcs("%s_0.1.trim_Cout" % TAMO_file)
    print "topm",topm
    
    ##
    # Writes the top motif to its own directory.
    if topm!="None":

        newdic={}
        ml=MotifTools.txt2motifs("%s_0.1.trim" % TAMO_file)

        for m in ml:
            
            if m.oneletter == topm:
                newdic[m.oneletter] = m

        save_motifs(newdic.values(),"%s.TOP" % TAMO_file)
        os.system("rm %s_0.1.trim" % TAMO_file)
        os.system("rm %s_0.1.trim_Cout" % TAMO_file)
    
    ##
    # Writes a blank document if there was no top motif.
    else:
        oup=open("%s.TOP" % TAMO_file,"w")
        oup.close()

def run_UPGMA(TAMO_file, height, wdir, walltime = 120, mem = '124'):
    
    '''This script creates a tree for motifs based on their PCC distance.
    
    The function works on the output of the combine_distance_matrix and 
    check_square_matrix function. The tree that is output is used by the 
    merge_runs function.
    
    The construction of the tree is performed by an R script, and it is submited
    as a job.
    
    This function is based on a script writen by Cheng Zou, and then converted
    to a function by Alex Seddon.
    '''
    
    # Get the path to pcc_merge_CC.py directory
    script_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
    
    w_hour = str(walltime / 60)
    w_min  = str(walltime % 60)
    
    job = open("%s_%s.UPGMA.sh" % (TAMO_file, height), 'w')
    job.write("#!/bin/sh --login\n\n")
    job.write("#PBS -q main\n")
    job.write("#PBS -l nodes=1:ppn=1,walltime=%s:%s:00,mem=%sgb\n" % (w_hour, w_min, mem))
    
    job.write("R --vanilla --slave --args %s/%s.dm  %s < %s/UPGMA_final.R > %s.Rout" % (wdir,TAMO_file,height,script_dir,TAMO_file))
    
    job.close()
    
    os.system("qsub %s_%s.UPGMA.sh" % (TAMO_file, height))
    

def merge_runs(TAMO_file, wdir, height, distance, ancestor, target, genome):
    
    '''This script is used to merge motifs with the PCC matrix of all motifs.
    
    The script was originally written by Cheng Zou, and then converted to a 
    function by Alex Seddon.
    '''
    
    print "Here are the parameters you specified in this run "
    print "-tamo        %s" % TAMO_file
    print "-wdir        %s" % wdir
    print "-h        height to cut the tree, %s" % height
    print "-distance    %s" % distance
    print "-ancestor    %s" % ancestor
    print "-target    %s" % target
    print "-genome    %s" % genome
        
    if TAMO_file == '' or wdir == '':
        help()

    os.system("cd %s" % wdir)

    os.chdir(wdir)
    
    
    # This code was in the original clustering script. It has been taken out 
    # because the processes involved take too long and have been replaced by 
    # the matrix creation scripts and the run_UPGMA script.
    #if distance==0:
    #    os.system("python /mnt/home/seddonal/scripts/5_motif_merging/3.calculate_distance_matrix.py   -i %s --dfunc pccrange" % TAMO_file)
    #os.system("R --vanilla --slave --args %s.dm  %s< /mnt/home/seddonal/scripts/5_motif_merging/UPGMA_final.R> %s.Rout" % (TAMO_file,height,TAMO_file))

    cl_dic={}
    n=0
    
    # The file, TAMO_file.dm_UPGMA_Cl_0.05, is inorder of the motifs that appear
    # in the TAMO_file. If two motifs have the same number, they are considered 
    # a part of the same cluster.
    # This loop pulls the clustering information out of this file and creats
    # the dictionary cl_dic = {cluster_index:{motif_index:'1'}}
    for line in open("%s.dm_UPGMA_Cl_%s" % (TAMO_file,height),"r"):
        
        # Gets the clusterindex of this motif
        cl=line.strip()
        
        # Adds the cluster index if it has not been 
        if not cl_dic.has_key(cl):
            cl_dic[cl]={}
            
        cl_dic[cl][n]="1" # Adds the motif to that cluster
        n+=1              # Increases the motif index for the next motif

    #print cl_dic
    
    ml=MotifTools.txt2motifs(TAMO_file)
    old=[] # List of motifs that are the sole members of a cluster.
    
    # I think I can divide up this portion of the code to create a series 
    print ancestor,ancestor==0
    if ancestor==0:
    
        # This loop Looks at each cluster and attempts to merge the motifs
        # in the cluster if there are multiple motifs.
        for i in cl_dic.keys():
            
            print i,cl_dic[i]
            
            # If there are multiple motifs in the cluster, it merges the motifs
            if len(cl_dic[i])>1:
            
                # Adds all of the motifs in the cluster to an object called 
                # mlist.
                mlist=[] 
                for j in cl_dic[i]:
                    mlist.append(ml[j])
                    
                # Saves these motifs to there own TAMO file.
                save_motifs(mlist,"%s_sub_%s.tm" % (TAMO_file,i))
                
                # I am fairly certain that this process of converting to TF and
                # then returning it to TAMO format is only for keeping the names 
                # consistent. I need to verify this suspicion
                tamo2tf("%s_sub_%s.tm" % (TAMO_file,i))
                os.system("cat  %s_sub_%s.tm.tf > %s_sub_%s_sum.tm.tf" % (TAMO_file,i,TAMO_file,i))
                tf2tamo("%s_sub_%s_sum.tm.tf" % (TAMO_file,i))

                # Gets the top motif in the cluster.
                pick_chunk_score(wdir, '%s_sub_%s_sum.tm.tf.tm' % (TAMO_file,i), target, genome)
                
                # Removes the files that were created.
                os.system("rm  %s_sub_%s_sum.tm.tf.tm" % (TAMO_file,i))
                os.system("rm %s_sub_%s_sum.tm.tf" % (TAMO_file,i))
                os.system("rm -R %s_sub_%s.tm.tf_ST*" % (TAMO_file,i))
            
            # If there is only one motif in the cluster, it leaves it alone, 
            # And adds it to old
            else:
                key=cl_dic[i].keys()[0]
                old.append(ml[key])
            
    

    if ancestor==1:
    
        # This loop Looks at each cluster and attempts to merge the motifs
        # in the cluster if there are multiple motifs.    
        for i in cl_dic.keys():
        
            print i,cl_dic[i]

            # If there are multiple motifs in the cluster, it merges the motifs            
            if len(cl_dic[i])>1:

                # Adds all of the motifs in the cluster to an object called 
                # mlist.
                mlist=[]
                for j in cl_dic[i]:
                    mlist.append(ml[j])
                
                # Saves these motifs to there own TAMO file.                
                save_motifs(mlist,"%s_sub_%s.tm" % (TAMO_file,i))
                
                # Merges the motifs in the same cluster using STAMP
                tamo2tf("%s_sub_%s.tm" % (TAMO_file,i))
                
                # Gets the JASPER motifs that best match the motifs from within
                # the cluster.
                os.system("STAMP -tf  %s_sub_%s.tm.tf  -sd /home/chengzou/bin/STAMP/ScoreDists/JaspRand_PCC_SWU.scores  \
                 -go  1000 -ge 1000 -cc PCC -align SWU -out %s_sub_%s.tm.tf_STout -chp > %s_sub_%s.tm.tf_STout.log" % (TAMO_file,i,TAMO_file,i,TAMO_file,i))
                parse_out_STAMP(TAMO_file,i)
                
                # combines the JASPER motifs with the cluster motif and then
                # converts them all to one TAMO file
                os.system("cat  %s_sub_%s.tm.tf %s_sub_%s.tm.tf_SToutFBP.txt.mod %s_sub_%s.tm.tf_STout_tree_clusters.txt > %s_sub_%s_sum.tm.tf" % (TAMO_file,i,TAMO_file,i,TAMO_file,i,TAMO_file,i))
                tf2tamo("%s_sub_%s_sum.tm.tf" % (TAMO_file,i))

                # Gets the top motif within the TAMO file.
                pick_chunk_score(wdir, '%s_sub_%s_sum.tm.tf.tm' % (TAMO_file,i), target, genome)
                
                # Removes any files created in the processing.
                os.system("rm  %s_sub_%s_sum.tm.tf.tm" % (TAMO_file,i))
                os.system("rm %s_sub_%s_sum.tm.tf" % (TAMO_file,i))
                os.system("rm -R %s_sub_%s.tm.tf_ST*" % (TAMO_file,i))
            else:
                key=cl_dic[i].keys()[0]
                old.append(ml[key])

    # Combine together the top motifs from every 
    os.system("cat %s_sub_*_sum.tm.tf.tm.TOP > %s_sub_new.tm" % (TAMO_file,TAMO_file))
    save_motifs(old,"%s_sub_old.tm" % (TAMO_file))
    os.system("cat %s_sub_old.tm %s_sub_new.tm > %s_P1.tm" % (TAMO_file,TAMO_file,TAMO_file))


def merge_runs_cc(TAMO_file, wdir, height, distance, ancestor, target, genome):
    
    '''This script is used to merge motifs with the PCC matrix of all motifs.
    
    The script was originally written by Cheng Zou, and then converted to a 
    function by Alex Seddon.
    '''
    
    print "Here are the parameters you specified in this run "
    print "-tamo        %s" % TAMO_file
    print "-wdir        %s" % wdir
    print "-h        height to cut the tree, %s" % height
    print "-ancestor    %s" % ancestor
    print "-target    %s" % target
    print "-genome    %s" % genome
        
    if TAMO_file == '' or wdir == '':
        help()

    os.system("cd %s" % wdir)

    os.chdir(wdir)
    
    # Get the directory where the script is located.
    script_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
    
    
    # This code was in the original clustering script. It has been taken out 
    # because the processes involved take too long and have been taken up by 
    # the matrrix creation scripts and the run_UPGMA script.
    #if distance==0:
    #    os.system("python /mnt/home/seddonal/gil scottscripts/5_motif_merging/3.calculate_distance_matrix.py   -i %s --dfunc pccrange" % TAMO_file)
    #os.system("R --vanilla --slave --args %s.dm  %s< /mnt/home/seddonal/scripts/5_motif_merging/UPGMA_final.R> %s.Rout" % (TAMO_file,height,TAMO_file))

    cl_dic={}
    n=0
    
    # The file, TAMO_file.dm_UPGMA_Cl_0.05, is inorder of the motifs that appear
    # in the TAMO_file. If two motifs have the same number, they are considered 
    # a part of the same cluster.
    # This loop pulls the clustering information out of this file and creats
    # the dictionary cl_dic = {cluster_index:{motif_index:'1'}}
    for line in open("%s.dm_UPGMA_Cl_%s" % (TAMO_file,height),"r"):
        
        # Gets the clusterindex of this motif
        cl=line.strip()
        
        # Adds the cluster index if it has not been 
        if not cl_dic.has_key(cl):
            cl_dic[cl]={}
            
        cl_dic[cl][n]="1" # Adds the motif to that cluster
        n+=1 # Increases the motif index for the next motif

    #print cl_dic
    
    ml=MotifTools.txt2motifs(TAMO_file)
    old=[] # List of motifs that are the sole members of a cluster.
    
    # I think I can divide up this portion of the code to create a series 
    print ancestor,ancestor==0
    
    cc_output = open('merge_runs_cc', 'w')
    
    if ancestor==0:
    
        # This loop Looks at each cluster and attempts to merge the motifs
        # in the cluster if there are multiple motifs.
        for i in cl_dic.keys():
            
            print i,cl_dic[i]
            
            # If there are multiple motifs in the cluster, it merges the motifs
            if len(cl_dic[i])>1:
            
                # Adds all of the motifs in the cluster to an object called 
                # mlist.
                mlist=[] 
                for j in cl_dic[i]:
                    mlist.append(ml[j])
                    
                # Saves these motifs to there own TAMO file.
                save_motifs(mlist,"%s_sub_%s.tm" % (TAMO_file,i))
                
                cc_output.write('module load TAMO; python %s/pcc_merge_CC.py merge_runs_no_ancestor -t %s/%s -i %s -target %s -genome %s\n' % (script_dir, wdir, TAMO_file, i, target, genome))
                
            # If there is only one motif in the cluster, it leaves it alone, 
            # And adds it to old
            else:
                key=cl_dic[i].keys()[0]
                old.append(ml[key])

    if ancestor==1:
    
        # This loop Looks at each cluster and attempts to merge the motifs
        # in the cluster if there are multiple motifs.    
        for i in cl_dic.keys():
        
            print i,cl_dic[i]

            # If there are multiple motifs in the cluster, it merges the motifs            
            if len(cl_dic[i])>1:

                # Adds all of the motifs in the cluster to an object called 
                # mlist.
                mlist=[]
                for j in cl_dic[i]:
                    mlist.append(ml[j])
                
                # Saves these motifs to there own TAMO file.                
                save_motifs(mlist,"%s_sub_%s.tm" % (TAMO_file,i))
                
                cc_output.write('module load TAMO; module load STAMPmotif; python %s/pcc_merge_CC.py merge_runs_ancestor -t %s/%s -i %s -target %s -genome %s\n' % (script_dir, wdir, TAMO_file, i, target, genome))
                
            else:
                key=cl_dic[i].keys()[0]
                old.append(ml[key])

    # Combine together the motifs that are in there own cluster.
    #os.system("cat %s_sub_*_sum.tm.tf.tm.TOP > %s_sub_new.tm" % (TAMO_file,TAMO_file))
    save_motifs(old,"%s_sub_old.tm" % (TAMO_file))
    #os.system("cat %s_sub_old.tm %s_sub_new.tm > %s_P1.tm" % (TAMO_file,TAMO_file,TAMO_file))

def merge_runs_no_ancestor(wdir, TAMO_file, i, target, genome):
    
    # I am fairly certain that this process of converting to TF and
    # then returning it to TAMO format is only for keeping the names 
    # consistent. I need to verify this suspicion, AS.
    tamo2tf(wdir, "%s/%s_sub_%s.tm" % (wdir, TAMO_file, i))
    os.system("cat  %s/%s_sub_%s.tm.tf > %s/%s_sub_%s_sum.tm.tf" % (wdir,TAMO_file,i,wdir,TAMO_file,i))
    tf2tamo(wdir, "%s_sub_%s_sum.tm.tf" % (TAMO_file,i))

    # Gets the top motif in the cluster.
    pick_chunk_score(wdir, '%s/%s_sub_%s_sum.tm.tf.tm' % (wdir,TAMO_file,i), target, genome)
    
    # Removes the files that were created.
    os.system("rm %s_sub_%s_sum.tm.tf.tm" % (wdir,TAMO_file,i))
    os.system("rm %s/%s_sub_%s_sum.tm.tf" % (wdir,TAMO_file,i))
    os.system("rm -R %s/%s_sub_%s.tm.tf_ST*" % (wdir,TAMO_file,i))


def merge_runs_ancestor(wdir, TAMO_file, i, target, genome):
    
    '''Merges the motifs within a single cluster.
    
    This function will identify motifs that are within the JASPER data set that
    are similar to the motifs within the cluster.
    '''

    # Merges the motifs in the same cluster using STAMP
    print "%s_sub_%s.tm" % (TAMO_file, i)
    tamo2tf("%s_sub_%s.tm" % (TAMO_file, i))
    
    # Gets the JASPER motifs that best match the motifs from within
    # the cluster.
    print "STAMP -tf  %s_sub_%s.tm.tf  -sd /home/chengzou/bin/STAMP/ScoreDists/JaspRand_PCC_SWU.scores -go  1000 -ge 1000 -cc PCC -align SWU -out %s_sub_%s.tm.tf_STout -chp > %s_sub_%s.tm.tf_STout.log" % (TAMO_file,i,TAMO_file,i,TAMO_file,i)
    os.system("STAMP -tf  %s_sub_%s.tm.tf  -sd /home/chengzou/bin/STAMP/ScoreDists/JaspRand_PCC_SWU.scores  \
     -go  1000 -ge 1000 -cc PCC -align SWU -out %s_sub_%s.tm.tf_STout -chp > %s_sub_%s.tm.tf_STout.log" % (TAMO_file,i,TAMO_file,i,TAMO_file,i))
    parse_out_STAMP_job(TAMO_file, i)
    
    # combines the JASPER motifs with the cluster motif and then
    # converts them all to one TAMO file
    os.system("cat  %s_sub_%s.tm.tf %s_sub_%s.tm.tf_SToutFBP.txt.mod %s_sub_%s.tm.tf_STout_tree_clusters.txt > %s_sub_%s_sum.tm.tf" % (TAMO_file,i,TAMO_file,i,TAMO_file,i,TAMO_file,i))
    tf2tamo("%s_sub_%s_sum.tm.tf" % (TAMO_file,i))

    # Gets the top motif within the TAMO file.
    pick_chunk_score(wdir, '%s_sub_%s_sum.tm.tf.tm' % (TAMO_file,i), target, genome)
    
    # Removes any files created in the processing.
    os.system("rm %s_sub_%s_sum.tm.tf.tm" % (TAMO_file,i))
    os.system("rm %s_sub_%s_sum.tm.tf"    % (TAMO_file,i))
    os.system("rm -R %s_sub_%s.tm.tf_ST*" % (TAMO_file,i))

def check_top_files(wdir):

    '''Check that each cluster has a TOP file.'''

    # Gets all of the files in the working directory.
    wdir_file_list = os.listdir(wdir)
    
    # Find which clusters that have a command line in the merge_runs_cc and 
    # determines which job file this cluster coresponds to.
    merge_runs_cc_file = open('%s/merge_runs_cc' % wdir, 'r')
    
    # Dictionary contianing the cluster information
    # cluster_cc_dict = {cluster number: (job number, command line)}
    cluster_cc_dict = {}
    job_number = 1
    for line in merge_runs_cc_file:
        cluster = line.strip().split('-i')[1].split(' ')[1]
        cluster_cc_dict[cluster] = (str(job_number), line)
        job_number += 1
    
    
    # Identify which clusters have a TOP file.
    
    cluster_TOP_dict = {} # Dictionary containing clusters with TOP files
    
    for file_name in wdir_file_list:
        if file_name.endswith('.TOP'):
            
            cluster = file_name.strip().split('_sub_')[1].split('_')[0]
            
            cluster_TOP_dict[cluster] = 1
            
    
    print 'There are', len(cluster_cc_dict), 'clusters in merge_runs_cc,'
    print 'and', len(cluster_TOP_dict), 'of the clusters have a TOP file.'
        
    output = open('failed_merge_runs_cc', 'w')
    for cluster in cluster_cc_dict:
        if not cluster in cluster_TOP_dict:
            print 'job'+cluster_cc_dict[cluster][0]+'.sh'
            output.write(cluster_cc_dict[cluster][1])            
    output.close()
    
            

def merge_runs_cleanup(wdir, TAMO_file):

    '''Concatonates all the motif files and deletes unneeded files.'''
    
    os.system("cd %s" % wdir)
    os.chdir(wdir)
    
    # Combine all of the TOP motif files and the single cluster motifs
    os.system("cat %s/%s_sub_*_sum.tm.tf.tm.TOP > %s/%s_sub_new.tm" % (wdir,TAMO_file,wdir,TAMO_file))
    os.system("cat %s/%s_sub_old.tm %s/%s_sub_new.tm > %s/%s_P1.tm" % (wdir,TAMO_file,wdir,TAMO_file,wdir,TAMO_file))

    # Delete the excess files. This is meant to be conservative so that you do
    # not lose the TOP files until you are sure you want to delete them. 
    os.system("rm %s/*.tf" % wdir)
    os.system('rm %s/*_SToutFBP.txt' % wdir)
    os.system('rm %s/*.mod' % wdir)
    os.system('rm %s/*_STout.log' % wdir)
    os.system('rm %s/*.tree' % wdir)
    os.system('rm %s/*_STout_tree_clusters.txt' % wdir)
    os.system('rm %s/*tm.tf.tm' % wdir)
    os.system('rm %s/*_0.1.trim' % wdir)
    os.system('rm %s/*_Cout' % wdir)

# help functions
def help():

    print 'python pcc_merge_CC.py [Function] -w working directory -t TAMO_file\n'
     
    print '~Functions'
    print '\tcreate_cc: Splits the TAMO file, and creates command lines for'
    print '\tgenerating PCC score matricies. Requires: -w, -t\n'
    print '\tcheck_jobs: You can run this script after you have run all of the' 
    print '\tPCC matrix creation jobs. It creates a failed_jobs_cc, which you can' 
    print '\trerun by using qsub_hpc.py. The funciton will print the errors for'
    print '\tthe failed jobs. Requires: -w\n'
    print '\tcombine_distance_matrix: After all jobs have run successfully, you can'
    print '\tuse this function to combine all of the distance matricies that were'
    print '\tcreated. Requires: -t, -w\n'
    print '\tcheck_square_matrix: Once the combined matrix is created, you can' 
    print '\tuse this function to check the matrix was formed correctly.'
    print '\tRequires: -d, -w\n'
    print '\trun_UPGMA: This script submits a job that creates the UPGMA tree'
    print '\tfrom the PCC score matrix. You need to run this script before'
    print '\trunning the merge_runs script. Requires: -t, -w. Options: -h,'
    print '\t-wall, -m\n'
    print '\tmerge_runs: Once the PCC matrix is combined and checked, this function'
    print '\twill create a tree for the motifs, cut the tree at the spcified height'
    print '\tand then merge the motifs bellow this height. Requires: -t, -w.'
    print '\tOptions: -h, -dist, -ancestor, -target, -genome\n'
    print '\tmerge_runs_cc: Divides the merge runs process into a series of jobs'
    print '\tthat you can run using qsub_hpc.py. Requires the same input as'
    print '\tmerge_runs\n'
    print '\tmerge_runs_ancestor and merge_runs_no_ancestor: Used to run individual'
    print '\tjobs from merge_runs_cc.\n'
    print '\tcheck_top_files: This function checks if each cluster that was'
    print '\tsupposed to be merged has a .TOP output file. The script reports'
    print '\twhich jobs the cluster belongs to. Requires: -w\n'
    print '\tmerge_runs_cleanup: After running merge_runs jobs, you can use this'
    print '\tto combine all of the TOP files and the single cluster motifs. It also'
    print '\tremoves all of the unnecessary files from generated in the merging'
    print '\tprocess. Requires: -t, -w\n'
    
    print '~Parameters'

    print "\t-ancestor := consider estimated ancestors of motifs(1) or not(0),"
    print '\t\tdefault is 1.'
    print '\t-d := Merged PCC matrix'
    print "\t-genome := As default, it is fasta fiand determines which job file"
    print "\t\ton the array and with promoter seqences annotated in TAIR8"
    print "\t-h := height to cut the tree, Be aware of when the height is equal"
    print "\t\tto 0.05,the distance between two taxon is 0.10"
    print "\t-mem  := memory required to run UPGMA. Default is 256 gb."
    print '\t-t := TAMO file that you wish to merge'
    print "\t-target := the list of target genes that you are interested,"
    print '\t\tas default it is all stress responsive genes q<0.05,2 fold'
    print '\t-w := Working directory for the merging process'
    print "\t-wall := walltime for running UPGMA. Default is 124 minutes."
    print '\t-c := Command line file'    

# Main
def main(): 
        
    # Get the function
    try:
        function = sys.argv[1]
    except IndexError:
        print __doc__
        sys.exit()
    
    # Set defaults
    TAMO_file = ''
    wdir      = os.path.abspath(os.curdir)
    distance  = 1
    ancestor  = 1
    height    = 0.05
    genome    = "/home/chengzou/project/Motif_run/_tamo_rawdata/all_on_array_withpromoter_seq.fa"
    target    = "/home/chengzou/project/Motif_run/_tamo_rawdata/expr_any_condition_withpromoter_seq"
    mem       = '124'
    walltime  = 120
    write_failed_cc = False
    dfunc = 'pccrange'
    
    # Get arguments
    for n in range(2, len(sys.argv), 2):
        if sys.argv[n] == '-w':
            wdir        = sys.argv[n+1].rstrip('/') 
        elif sys.argv[n] == '-t':
            TAMO_file     = sys.argv[n+1]
        elif sys.argv[n] == '-t2':
            TAMO_file_2   = sys.argv[n+1]
        elif sys.argv[n] == '-d':
            PCC_matrix    = sys.argv[n+1]
        elif sys.argv[n] == "-h":
            height        = float(sys.argv[n+1])
        elif sys.argv[n] == "-dist":
            distance      = int(sys.argv[n+1])     
        elif sys.argv[n] == "-ancestor":
            ancestor      = int(sys.argv[n+1])
        elif sys.argv[n] == "-target":
            target        = sys.argv[n+1]
        elif sys.argv[n] == '-m':
            mem           = sys.argv[n+1]
        elif sys.argv[n] == '-wall':
            walltime      = int(sys.argv[n+1])
        elif sys.argv[n] == '-i':
            i             = sys.argv[n+1]
        elif sys.argv[n] == '-genome':
            genome        = sys.argv[n+1]
        elif sys.argv[n] == '-c':
            cc_file       = sys.argv[n+1]
        elif sys.argv[n] == '--dfunc':
            dfunc = sys.argv[n+1]
            
    
    if function == 'create_cc':
        create_cc(wdir, TAMO_file, dfunc)
    elif function == 'create_cc_2':
        create_cc_for_two(wdir, TAMO_file, TAMO_file_2)
    elif function == 'check_jobs':
        check_jobs(wdir)
    elif function == 'check_outputs':
        check_outputs(wdir, cc_file)
    elif function == 'combine_distance_matrix':
        combine_distance_matrix(wdir, TAMO_file)
    elif function == 'check_square_matrix':
        check_square_matrix(wdir, PCC_matrix)
    elif function == 'combine_distance_matrix_2':
        combine_distance_matrix_for_2(wdir, TAMO_file, TAMO_file_2)
    elif function == 'run_UPGMA':
        run_UPGMA(TAMO_file, height, wdir, walltime, mem)
    elif function == 'merge_runs':
        merge_runs(TAMO_file, wdir, height, distance, ancestor, target, genome)
    elif function == 'merge_runs_cc':
        merge_runs_cc(TAMO_file, wdir, height, distance, ancestor, target, genome)
    elif function == 'merge_runs_no_ancestor':
        merge_runs_no_ancestor(wdir, TAMO_file, i, target, genome)
    elif function == 'merge_runs_ancestor':
        merge_runs_ancestor(wdir, TAMO_file, i, target, genome)
    elif function == 'merge_runs_cleanup':
        merge_runs_cleanup(wdir, TAMO_file)
    elif function == 'check_top_files':
        check_top_files(wdir)
    elif function == 'help':
        print '/'.join(os.path.abspath(__file__).split('/')[:-1])
        print __doc__
        sys.exit()
    else:
        print __doc__
        sys.exit()

if __name__ == '__main__':
    main()

