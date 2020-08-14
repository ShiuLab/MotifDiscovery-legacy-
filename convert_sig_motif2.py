##script used to parse PCC motif matrix to sig or not sig
import os, sys

sum_matrix = open(sys.argv[1],"r") #miteupLMup_TFBM_TFfamily_Low_PCC_average.txt
output = open(sys.argv[1] + "_sig.txt", "w")


D = {}
def add_data_to_dict(inp):
    for line in inp:
        if line.startswith("#python"):
            pass
        if line.startswith("motif"):
            L = line.strip().split("\t")
            title_list = L[0:]
            title_str = '\t'.join(title_list)
            output.write(title_str + '\n')
        else:
            L2 = line.strip().split('\t')
            motif = L2[0]
            PCC_dist = L2[1:]
            sig_list= []
            for dist in PCC_dist:
                if float(dist) <=0.39:
                    sig_list.append(str(1))
                else:
                    sig_list.append(str(0))
            sig_str = "\t".join(sig_list)
            output.write('%s\t%s\n' % (motif, sig_str))

add_data_to_dict(sum_matrix)

output.close()