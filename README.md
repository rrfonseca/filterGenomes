# from assemblies to Reference Genomes

###sums length of fasta seq
###assumes no missing seqs and no missing lines
###sort sca by size, largest first
###rename sca accordingly

fa=assembly.fasta
tag=fasta
prefix=${fa%$tag}

python sort-assembly-rename-sca.py $fa 

###index genome

./index-genome.sh $fa

###get assembly stats with https://github.com/WenchaoLin/assemblyStatics

./assemblyStatics.py $fa > ${fa}.stats

###get bed file with chr size to create the non-repeat coor file

for i in *oldIds.sizes ; do awk '{print $2}' $i > ${i}.a ; done
for i in *oldIds.sizes ; do awk '{print 0}' $i > ${i}.b ; done
for i in *oldIds.sizes ; do awk '{print $3}' $i > ${i}.c ; done
for i in *oldIds.sizes ; do paste ${i}.a ${i}.b ${i}.c > ${i}.chr.bed ; rm ${i}.a ${i}.b ${i}.c ; done

###calculate mappability with GENMAP https://github.com/cpockrandt/genmap

##create the index, only has to be done once

$genmap index -F ${fa}.simpleIds.sortSize.fa -I ${prefix}.genmap.out >& ${prefix}.genmap.log

##calculate mappability and export bedgraph file
$genmap map -K 100 -E 2 -I ${prefix}.genmap.out -O ${prefix}.mappability -bg -T 8 >& ${prefix}.genmap.map.log

###detect repeat regions with RepeatMasker https://www.repeatmasker.org/

$RM -pa 40  -frag 50000 -species chicken ${fa}.simpleIds.sortSize.fa -dir ./ -xm -xsmall -gff >& ${prefix}.repeatMasker.chicken.log

###detect potential sex-linked scaffolds with SATC https://github.com/popgenDK/SATC
###requires list of bamfiles

Rscript --vanilla $satc -i bamList -o ${prefix}.median.m10k.M10 --useMedian TRUE --minLength 10000 --M 10

###intersect output from SATC XZ_scaff.list that contains putative sex-linked or weird depth scaffolds
###with the bed file containing the full length scaffold coors

awk 'FNR==NR{A[$1]=1;next} A[$1]' ${prefix}.median.m10k.M10_XZ_scaff.list ${fa}.oldIds.sizes.chr.bed > ${prefix}.m10k.M10_XZ_scaff.bed

###get GENMAP areas of mappability =1 (*bedgraph)

awk '$4==1 {print $0}' ${prefix}.mappability.bedgraph | cut -f1-3 > ${prefix}.mappability.1.bedgraph

###next steps require bedtools https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html

###remove repeat intervals (RepeatMasker output *out.gff )

bedtools subtract -a ${prefix}.mappability.1.bedgraph -b ${prefix}.fa.simpleIds.sortSize.fa.out.gff > ${prefix}.map1.notRep.bed

###get bad sites for masking fasta by subtracting good bed from full bed

bedtools subtract -a ${prefix}.fa.oldIds.sizes.chr.bed -b ${prefix}.map1.notRep.bed > ${prefix}.notMap1.withRep.bed

###mask fasta genome

bedtools maskfasta -fi $fa -fo ${prefix}.masked.fa -bed ${prefix}.notMap1.withRep.bed

###index-genome.sh

index-genome.sh ${prefix}.masked.fa  

