#sums length of fasta seq
#assumes no missing lines
#calculates p missing (Ns)

import sys
import re

inFasta=sys.argv[1]

out=open("%s.sizes"%(inFasta),"w")

infos=open(inFasta)
info=infos.readline()
n=0
miss=0
while info:
	if re.match(">",info):
		if n>0:
			out.write("%s\t%s\n"%(n,miss*100/n*0.01))
		out.write("%s\t"%(info.rstrip()[1:]))
		miss=0
		n=0
	else:
		n+=len(info.rstrip())
		miss+=info.count("N")

	info=infos.readline()

out.write("%s\t%s\n"%(n,miss*100/n*0.01))

out.close()
