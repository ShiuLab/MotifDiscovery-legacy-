#calculate the average PCC distance within the TF families; outpu the lowest
#example: python /Users/liumingjung/python_script/TF_with_between_correlation_low_PCC_average.py 
#input [1] = ../Weirauch_2014/Athaliana_TFBM_v1.01.tm.index.direct.index_subset.txt 
#input [2] =../WFup.txt.tamo-Athaliana_TFBM_v1.01.tm.index.direct.index.tm.dm

import sys

#input the order of TFs within/between families
#/Users/liumingjung/Desktop/Ara_coronatine_20130905/pCRE_analyses/correaltion_know_motif/Weirauch_2014/Athaliana_TFBM_v1.01.tm.index.direct.index_subset.txt
in1=open(sys.argv[1], 'r')
tf={}
motif=[]
for l1 in in1:
    if l1.startswith("#"):
          #line=line.strip("\n").strip("\r").strip()
          #out.write("%s\n"%(line))
          #print "test", line
        pass
    elif l1=="\n":
        pass
    else:

        l=l1.strip("\n").strip("\r").strip()
        p=l.strip().split("\t")
        tf[p[0]+"_"+p[3]]=p[3]
        motif.append(p[3])

print "TF motif number:", len(motif), motif
in1.close()
 

out=open(sys.argv[2]+"_TFs_Low_PCC_average.txt", 'w')
out.write('#python %s\n'%' '.join(sys.argv))
for i in set(tf.values()):
    out.write("%s\t"%(i))
out.write("\n")


#input the PCC distance metrix for each comparison
in3=open(sys.argv[2], 'r')


for l3 in in3:
    if l3.startswith("#"):
        pass
    elif l3=="\n":
        pass
    else:
              #lo=[]
        ll=l3.strip("\n").strip("\r").strip()
        pp=ll.strip().split("\t")
        print "comparison motif number:", len(pp)
        motif_index={}
        for i in range(0, len(pp)):
            s=str(i)+"_"+pp[i]
            motif_index[s]=motif[i]
        for m in set(tf.values()):
            lo=[]
            sum=0
            n=0
            print m
            for i, j in motif_index.iteritems():
                if j ==m:
                    print j
                    p=i.split("_")
                    sum=sum+float(p[1])
                    lo.append(float(p[1]))
                    n+=1
            av=sum/n
            out.write("%s\t"%(av))
            print m, lo
    out.write("\n")
print "number_input_motif:",len(motif)

in3.close()
out.close()
