#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=4
#SBATCH --output=./fastp.out
#SBATCH -e ./fastp.err

source /storage/herman/home/xwen/miniconda3/etc/profile.d/conda.sh

#This code assumes you are in some working directory of arbitrary name.
#You should have an illumina run folder in that directory. With a sample sheet present.


WD=/storage/herman/home/xwen/RNAseq/Runs
Run=/practice_fastp

conda activate fastp

cd $WD$Run

for i in $(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename "$F" | rev | cut -c 13- | rev; done | sort | uniq); do
    echo "Processing file: $i"
    srun \
	--nodes=1 \
	--mem=10GB \
	--ntasks=1 \
	--cpus-per-task=4 \
	fastp --in1 "$i"_R1.fastq.gz --in2 "$i"_R2.fastq.gz --out1 "$i"_R1_trim.fastq.gz --out2 "$i"_R2_trim.fastq.gz --thread 4 --trim_poly_g --poly_g_min_len=11 --length_required 25 --html $WD$Run/2_FASTQ/"$i".fastp.html --json $WD$Run/2_FASTQ/"$i".fastp.jzon
done

echo "finish"
