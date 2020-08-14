import os,sys
import pandas as pd

dir1= sys.argv[1] # input the starting directory with TFs_Low_PCC_average.txt files or distance files
type1= str(sys.argv[2]) #1= TFs_Low_PCC_average files, 2= distance file (.dm)

def combine_files(motif, df1):
    #df1['Motif'] = motif
    print("motifs: ", len(motif))
    print("rows in df: ", len(df1.index))
    if df1.empty == True:
        print("df1 is empty")
    else:
        df1.insert(0, "Motif", motif, True)
    return df1

def get_TF_header(inp):
    TF_list= []
    for line in inp:
        L=line.strip().split("\t")
        TF_list.append(str(L[0]))
    return (TF_list)

if type1 == "1":
    for file in os.listdir(dir1):
        if file.endswith("_TFs_Low_PCC_average.txt"):
            name = file.strip().split(".tamo")[0]
            print (name)
            df1 = pd.read_csv((dir1 + "/" + file), sep='\t', header=1)
            inp = open(dir1 + "/" + name)
            motif_list=[]
            for line in inp:
                L=line.strip().split()
                motif_list.append(L[0])
            newdf= combine_files(motif_list, df1)
            newdf.to_csv(path_or_buf=str(file)+"_motifs.txt", sep="\t", header=True)
            inp.close()
 
elif type1 == "2":
    TF_file = open(sys.argv[3], 'r')
    TF_list= get_TF_header(TF_file)
    TF_file.close()
    for file in os.listdir(dir1):
        if file.endswith(".tm.dm"):
            name = file.strip().split(".tamo")[0]
            print (name)
            df1 = pd.read_csv((dir1 + "/" + file), sep='\t', header=None)   
            df1.columns = TF_list
            inp = open(dir1 + "/" + name)
            motif_list=[]
            for line in inp:
                L=line.strip().split()
                motif_list.append(L[0])
            newdf= combine_files(motif_list, df1)
            newdf.to_csv(path_or_buf=str(file)+"_motifs.txt", sep="\t", header=True)
            inp.close()
            
else:
    print ("need TYPE; 1= _TFs_Low_PCC_average.txt files, 2= .tm.dm. For 2, also need 3rd argument (TF index to add header)")
            