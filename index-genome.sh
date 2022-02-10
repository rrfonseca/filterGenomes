#index fasta genome file for bwa and GATK

fa=$1
tag=.fa
id=${fa%$tag}

bwa index -a bwtsw $1

picard-tools CreateSequenceDictionary REFERENCE=${fa} OUTPUT=${id}.dict

nice samtools faidx ${fa}


