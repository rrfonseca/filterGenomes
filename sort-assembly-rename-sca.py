#sums length of fasta seq
#assumes no missing seqs and no missing lines
#sort sca by size, largest first
#rename sca accordingly 

import sys
import re

inFasta=sys.argv[1]

out=open("%s.oldIds.sizes"%(inFasta),"w")
out1=open("%s.simpleIds.sortSize.fa"%(inFasta),"w")

infos=open(inFasta)
info=infos.readline()

seqD={}
sizeD={}

print "reading fasta\n" 

key=""
i=0
while info:
	if re.match(">",info):
		key=info.rstrip()[1:]
		seqD[key]=""
		sizeD[key]=0
		i+=1
	else:
		seqD[key]+=info
		sizeD[key]+=len(info.rstrip())
	print i		
	info=infos.readline()

print "done reading\n"

sizeL=[]
for key in sizeD:
	sizeL.append((sizeD[key],key))

sizeL.sort()

print "done sorting\n"

i=len(sizeL)-1
k=1
while i>-1:
	size=sizeL[i][0]
	key=sizeL[i][1]
	seq=seqD[key]
	out1.write(">s%s\n%s"%(k,seq))
	out.write("%s\ts%s\t%s\n"%(key,k,size))	
	k+=1
	i-=1

out.close()
out1.close()

