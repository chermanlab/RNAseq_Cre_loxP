#!/bin/bash

CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh

#This code assumes you are in some working directory of arbitrary name.
#You should have an illumina run folder in that directory. With a sample sheet present.

# mkdir 1_DEMULTIPLEXING 2_FASTQ 3_SALMON 4_DESEQ2
# conda activate bcl2fastq2
# mv ./*NB* ./0_RUNFOLDER
# cd ./0_RUNFOLDER
# bcl2fastq -R ./ -o ../ --no-lane-splitting --barcode-mismatches 1
# cd ../
# mv Reports Stats 1_DEMULTIPLEXING
# filenames_R1=($(find ./ -maxdepth 1 -type f -name "*R1*.fastq.gz" | while read F; do basename $F | grep -oP '.*(?=_S\d+)'; done | sort | uniq))
# filenames_R2=($(find ./ -maxdepth 1 -type f -name "*R2*.fastq.gz" | while read F; do basename $F | grep -oP '.*(?=_S\d+)'; done | sort | uniq))
# for i in "${filenames_R1[@]}"; do mv "$i"*R1*.fastq.gz "$i"_R1.fastq.gz; done
# for i in "${filenames_R2[@]}"; do mv "$i"*R2*.fastq.gz "$i"_R2.fastq.gz; done

# conda activate fastp
# max_execution_time=600 #50M 150BP SE Reads took 1200s, so 2400s gives adequate time
# for i in $(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename "$F" | rev | cut -c 13- | rev; done | sort | uniq); do
    # while true; do
        # echo "Processing file: $i"
        # setsid fastp --in1 "$i"_R1.fastq.gz --in2 "$i"_R2.fastq.gz --out1 "$i"_R1_trim.fastq.gz --out2 "$i"_R2_trim.fastq.gz --thread 12 --trim_poly_g --poly_g_min_len=11 --length_required 25 --html 2_FASTQ/"$i".fastp.html --json 2_FASTQ/"$i".fastp.jzon & fastp_pid=$!
        # (sleep $max_execution_time; pkill -P $fastp_pid 2>/dev/null; kill -9 $fastp_pid 2>/dev/null) &  # TIMER and Kill fastp after timeout
        # wait $fastp_pid # Wait for fastp to complete or for the timer to expire
        # fastp_exit_status=$?
        # if [ $fastp_exit_status -eq 0 ]; then # Check if fastp exited normally (0) or was killed by the timer (137)
			# mv "$i"*R?.fastq.gz ./2_FASTQ
            # break  # Exit the while loop and proceed to the next file
        # elif [ $fastp_exit_status -eq 137 ]; then
            # echo "fastp took too long or froze. Retrying file: $i"
        # else
            # echo "fastp encountered an error. Skipping file: $i"
            # break  # Exit the while loop and proceed to the next file
        # fi
    # done
# done

# conda activate salmon
# for i in $(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 18- | rev; done | sort | uniq)
	# do salmon quant \
	# --libType IU \
	# --index /mnt/e/TRANSCRIPTOMICS/SALMON_INDEXES/MG1655_CDS_K17 \
	# --mates1 "$i"_R1_trim.fastq.gz \
	# --mates2 "$i"_R2_trim.fastq.gz \
	# --output 3_SALMON/"$i"_salmon_CDS \
	# --seqBias \
	# --gcBias \
	# --threads 14 \
	# --allowDovetail
	# salmon quant \
	# --libType IU \
	# --index /mnt/e/TRANSCRIPTOMICS/SALMON_INDEXES/MG1655_ncRNA_K17 \
	# --mates1 "$i"_R1_trim.fastq.gz \
	# --mates2 "$i"_R2_trim.fastq.gz \
	# --output 3_SALMON/"$i"_salmon_ncRNA \
	# --seqBias \
	# --gcBias \
	# --threads 14 \
	# --allowDovetail
	# mv "$i"_R?_trim.fastq.gz ./2_FASTQ
# done;

cd ./3_SALMON
current_dir=$(pwd)
main_dir="$current_dir"
for dir in "$current_dir"/*_salmon_CDS/; do
    if [ -d "$dir" ]; then
        folder_name=$(basename "$dir" | sed 's/_salmon_CDS//')
        input_file="$dir/quant.sf"
        output_file="$current_dir/${folder_name}_quant_filtered.sf"
		awk '{print $1}' $input_file > Name.sf
        awk 'BEGIN {OFS="\t"} NR==1 {print "'"$folder_name"'"} NR>1 {print $4}' "$input_file" > "$output_file"
    fi
done
temp=$(find ./ -maxdepth 1 -type f -name "*.sf" | while read F; do basename "$F"; done)
temp1=$(echo -n "$temp" | tr '\n' ' ')
paste $temp1 > RNASEQ.sf 




