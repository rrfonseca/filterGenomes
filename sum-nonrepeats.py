#NOTE: sex sca will show up has having zero sites covered

# Xenodacnis_parina.map1.notRepMask.m10k.notXZ.bed
#s1      0       140
#s1      191     6194
#s1      6237    13912

# Xenodacnis_parina.fa.oldIds.sizes.chr.bed
#s1      0       811559
#s2      0       716767
#s3      0       561130

import sys

minScaSizeKb=int(sys.argv[1])

inOKSites=sys.argv[2] #"Xenodacnis_parina.map1.notRepMask.m10k.notXZ.bed"
inScaLength=sys.argv[3] #Xenodacnis_parina.fa.oldIds.sizes.chr.bed 

sp=inOKSites.split(".")[0]
outName="%s.map1.notRep.m%skb-OK.notXZ.sum.txt"%(sp,minScaSizeKb) #"Xenodacnis_parina.map1.notRepMask.m10k.notXZ.big%sk.sum.txt"%(minScaSizeKb)

out=open(outName,"w")

infos=open(inScaLength)
info=infos.readline()

scaL=[]
scaD={} #[sca]=finalSize

while info:
	x=info.rstrip().split()
	sca=x[0]
	size=int(x[2])
	scaD[sca]=0
	scaL.append((sca,size))

	info=infos.readline()

infos=open(inOKSites)
info=infos.readline()

while info:
	x=info.rstrip().split()
	sca=x[0]
	finalSize=scaD[sca]+(int(x[2])-int(x[1]))

	scaD[sca]=finalSize

	info=infos.readline()

out.write("scaffold\tsizeAssembly\tsizeNotRepeat\tpRepeats\n")
i=0
totalOK=0
total=0
while i<len(scaL):
	sca=scaL[i][0]
	fullSize=scaL[i][1]
	sizeOKsites=scaD[sca]

	pRepeats=round((fullSize-sizeOKsites)*100/fullSize*0.01,2)
	
	if sizeOKsites >= minScaSizeKb*1000:
		out.write("%s\t%s\t%s\t%s\n"%(sca,fullSize,sizeOKsites,pRepeats))
		totalOK+=sizeOKsites
	total+=fullSize

	i+=1

pRepeats=(total-totalOK)*100/total*0.01

print "Total size after filtering repeats, sex and <%sK: %sMb (fraction removed: %s)\n"%(minScaSizeKb,int(totalOK/1000000),pRepeats)

out.close()
