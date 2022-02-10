# prepareReferenceGenomesNonModel

#sums length of fasta seq
#assumes no missing seqs and no missing lines
#sort sca by size, largest first
#rename sca accordingly

fa=assembly.fasta
tag=fasta
sample=${fa%$tag}

python sort-assembly-rename-sca.py $fa 

#index genome

./index-genome.sh $fa

#get assembly stats with https://github.com/WenchaoLin/assemblyStatics

./assemblyStatics.py $fa > ${fa}.stats

#get bed file with chr size to create the non-repeat coor file

for i in *oldIds.sizes ; do awk '{print $2}' $i > ${i}.a ; done
for i in *oldIds.sizes ; do awk '{print 0}' $i > ${i}.b ; done
for i in *oldIds.sizes ; do awk '{print $3}' $i > ${i}.c ; done
for i in *oldIds.sizes ; do paste ${i}.a ${i}.b ${i}.c > ${i}.chr.bed ; rm ${i}.a ${i}.b ${i}.c ; done

###calculate mappability with GENMAP https://github.com/cpockrandt/genmap

#create the index, only has to be done once

$genmap index -F ${fa}.simpleIds.sortSize.fa -I ${sample}.genmap.out >& ${sample}.genmap.log

#calculate mappability and export bedgraph file
$genmap map -K 100 -E 2 -I ${sample}.genmap.out -O ${sample}.mappability -bg -T 8 >& ${sample}.genmap.map.log

###detect repeat regions with RepeatMasker https://www.repeatmasker.org/

$RM -pa 40  -frag 50000 -species chicken ${fa}.simpleIds.sortSize.fa -dir ./ -xm -xsmall -gff >& ${sample}.repeatMasker.chicken.log

###detect potential sex-linked scaffolds with SATC https://github.com/popgenDK/SATC
#requires list of bamfiles

Rscript --vanilla $satc -i bamList -o ${sample}.median.m10k.M10 --useMedian TRUE --minLength 10000 --M 10

###intersect output from SATC *XZ_scaff.list that contains putative sex-linked or weird depth scaffolds
#  with the bed file containing the full length scaffold coors

awk 'FNR==NR{A[$1]=1;next} A[$1]' ${sample}.median.m10k.M10_XZ_scaff.list ${fa}.oldIds.sizes.chr.bed > ${sample}.m10k.M10_XZ_scaff.bed

#get GENMAP areas of mappability =1 (*bedgraph)

awk '$4==1 {print $0}' ${sample}.mappability.bedgraph | cut -f1-3 > ${sample}.mappability.1.bedgraph

#next steps require bedtools https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html

#remove repeat intervals (RepeatMasker output *out.gff )
bedtools subtract -a ${sample}.mappability.1.bedgraph -b ${sample}.fa.simpleIds.sortSize.fa.out.gff > ${sample}.mappability1.notRepeatMasker.bed

#remove sex-linked sca from previous list
bedtools subtract -a ${sample}.mappability1.notRepeatMasker.bed -b ${sample}.m10k.M10_XZ_scaff.bed >${sample}.map1.notRepMask.m10k.notXZ.bed

#check final number on non-repeat/low map sites and output sca info if more than Xkb OK sites (non-repeat + map1 + not-sex)

python sum-nonrepeats.py 10 ${sample}.map1.notRepMask.m10k.notXZ.bed ${sample}.fa.oldIds.sizes

awk 'FNR==NR{A[$1]=1;next} A[$1]' ${sample}.map1.notRep.m10kb-OK.notXZ.sum.txt ${sample}.map1.notRepMask.m10k.notXZ.bed > ${sample}.map1.notRepMask.m10k-OKsites.notXZ.bed

#get list of sca names for PSMC

cut -f1 ${sample}.map1.notRepMask.m10k-OKsites.notXZ.bed | uniq > ${sample}.map1.notRepMask.m10k-OKsites.notXZ.sca

#mask fasta genome
#index masked genome





