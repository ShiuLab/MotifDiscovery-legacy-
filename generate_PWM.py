# generate PWM for indicated motifs
#module load TAMO (on hpc)
import sys
from TAMO import MotifTools
#input file: motif list
input=open(sys.argv[1], 'r')

motif=[]
for line in input:
    if line.startswith("#"):
        pass
    else:
        line=line.strip("\n").strip("\r").strip()
        motif.append(line)

input.close()

print (motif)

pw=[]
for i in range(0, len(motif)):
    m=MotifTools.Motif_from_text(motif[i])
    pw.append(m)
MotifTools.save_motifs(pw,sys.argv[1]+'.tamo')
