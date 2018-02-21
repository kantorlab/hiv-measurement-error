#!/bin/bash
set -e
i=${1:0:4}
j=${1:4:1}
CPU=$2

cd scratch/PrimerID/scripts/PreprocessedData

MIN_PHRED_SCORE=30

# 1) Perform full prinseq quality filtering
echo ${i}${j}
rm -rf ${i}${j}
mkdir ${i}${j}
cd ${i}${j}
echo -e "PRINSEQ filtering\n================="
prinseq-lite.pl \
	-fastq ../../RawData/${i}${j}R1.fastq \
	-fastq2 ../../RawData/${i}${j}R2.fastq \
	-out_format 3 \
	-out_good ${i}${j}_QC_1_R \
	-out_bad null \
	-min_len 230 \
	-trim_qual_left ${MIN_PHRED_SCORE} \
	-trim_qual_right ${MIN_PHRED_SCORE} \
	-trim_qual_window 5

# 2) Replace ALL low-quality nucleotides
echo -e "\n\nNucleotide masking\n=================="
for r in R_1 R_2
do
	fastq_masker \
		-q ${MIN_PHRED_SCORE} \
		-i ${i}${j}_QC_1_${r}.fastq \
		-o ${i}${j}_nucMask_2_${r}.fastq \
		-v
done
# Discard reads with missing mate
rm -rf *singletons*

cd ../../Alignments

echo ${i}${j}
rm -rf ${i}${j}
mkdir ${i}${j}
cd ${i}${j}

echo "Aligning after first QC step"

../../../pidalign \
	-t $CPU \
	-r1 ../../PreprocessedData/${i}${j}/${i}${j}_QC_1_R_1.fastq \
	-r2 ../../PreprocessedData/${i}${j}/${i}${j}_QC_1_R_2.fastq \
	-o ${i}${j}_QC_1.sam \
	--ref ../../References/5VM_${i}_N_withPrimer_cons.fasta \
	-v 1 \
	-s 210 \
	-d 4 \
	-Q

#echo "Aligning after second nucleotide masking step"
#
#../../../pidalign \
#	-t $CPU \
#	-r1 ../../PreprocessedData/${i}${j}/${i}${j}_nucMask_2_R_1.fastq \
#	-r2 ../../PreprocessedData/${i}${j}/${i}${j}_nucMask_2_R_2.fastq \
#	-o ${i}${j}_nucMask_2.sam \
#	--ref ../../References/5VM_${i}_N_withPrimer_cons.fasta \
#	-v 1 \
#	-s 210 \
#	-d 20

../../../pidalyse --r3223 ../../References/5VM_3223.fasta --r3236 ../../References/5VM_3236.fasta ${i}${j}_QC_1.sam

#../../../pidalyse --r3223 ../../References/5VM_3223.fasta --r3236 ../../References/5VM_3236.fasta ${i}${j}_nucMask_2.sam

cd "../../Comparing estimators"

echo ${i}${j}
rm -rf ${i}${j}
mkdir ${i}${j}
cd ${i}${j}

echo "Trim reads to fit and remove primerID simultaneously"
AmpliconClipper -i ../../Alignments/${i}${j}/${i}${j}_QC_1.sam -o ${i}${j}_clipped.sam -a ../${i}_amplicon.txt

echo "Extract reads with picard"
picard SamToFastq "I=${i}${j}_clipped.sam" "F=${i}${j}_R1.fastq" "F2=${i}${j}_R2.fastq" FU=single.fastq VALIDATION_STRINGENCY=LENIENT

