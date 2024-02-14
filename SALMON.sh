#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=4
#SBATCH --output=./salmon.out
#SBATCH -e ./salmon.err

source /storage/herman/home/xwen/miniconda3/etc/profile.d/conda.sh


WD=/storage/herman/home/xwen/RNAseq/Runs
Run=/practice_fastp

conda activate salmon

:<<COMMENT
cd $WD/$Run

mkdir "3_SALMON"

for i in $(find ./ -maxdepth 1 -type f -name "*trim.fastq.gz" | while read F; do basename $F | rev | cut -c 18- | rev; done | sort | uniq); do 
	mkdir "3_SALMON/"$i"_salmon_CDS" "3_SALMON/"$i"_salmon_NCRNA"
done

for i in $(find ./ -maxdepth 1 -type f -name "*trim.fastq.gz" | while read F; do basename $F | rev | cut -c 18- | rev; done | sort | uniq); do
	echo "Processing: $i"
	srun \
		--nodes=1 \
       		--mem=50GB \
        	--ntasks=1 \
        	--cpus-per-task=4 \
        	salmon quant \
        	--libType IU \
        	--index /storage/herman/home/xwen/RNAseq/INDEX_practice/Reference_genomes/_SALMON_INDEXES/MG1655/MG1655_CDS_K17 \
        	--mates1 "$i"_R1_trim.fastq.gz \
        	--mates2 "$i"_R2_trim.fastq.gz \
        	--output /storage/herman/home/xwen/RNAseq/Runs/practice_fastp/3_SALMON/"$i"_salmon_CDS \
        	--seqBias \
        	--gcBias \
        	--threads 14 \
        	--allowDovetail
	srun \
                --nodes=1 \
                --mem=50GB \
                --ntasks=1 \
                --cpus-per-task=4 \
        	salmon quant \
        	--libType IU \
        	--index /storage/herman/home/xwen/RNAseq/INDEX_practice/Reference_genomes/_SALMON_INDEXES/MG1655/MG1655_NCRNA_K17 \
        	--mates1 "$i"_R1_trim.fastq.gz \
        	--mates2 "$i"_R2_trim.fastq.gz \
       		--output /storage/herman/home/xwen/RNAseq/Runs/practice_fastp/3_SALMON/"$i"_salmon_NCRNA \
        	--seqBias \
        	--gcBias \
        	--threads 14 \
        	--allowDovetail 
       	mv "$i"_R?_trim.fastq.gz ./2_FASTQ
done
COMMENT

cd $WD/$Run/3_SALMON
current_dir=$(pwd)
main_dir="$current_dir"
for dir in "$current_dir"/*_salmon_CDS/; do
    if [ -d "$dir" ]; then
        folder_name=$(basename "$dir" | sed 's/_salmon_CDS//')
        input_file="$dir/quant.sf"
        output_file="$current_dir/${folder_name}_quant_filtered.sf"
                awk '{print $1}' $input_file > 0_Name.sf
        awk 'BEGIN {OFS="\t"} NR==1 {print "'"$folder_name"'"} NR>1 {print $4}' "$input_file" > "$output_file"
    fi
done
temp=$(find ./ -maxdepth 1 -type f -name "*.sf" | while read F; do basename "$F"; done | sort)
temp=( "${temp[@]/Undetermined_quant_filtered.sf/}" )
#temp1=$(echo -n "$temp" | tr '\n' ' ')
paste $temp > RNASEQ.sf


