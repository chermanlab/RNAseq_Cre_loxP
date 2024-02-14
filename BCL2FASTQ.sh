#!/usr/bin/sh
#SBATCH -t 3-00:00:00
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=4
#SBATCH --output=./bcl2cs.out
#SBATCH -e ./bcl2cs.err

source /storage/herman/home/xwen/miniconda3/etc/profile.d/conda.sh

conda activate RNAseq_1

cd /storage/herman/home/xwen/RNAseq/230915_11494a_11778a_pool/230915_NB551658_0272_AHWFJVBGXT/

bcl2fastq --help

srun bcl2fastq -R ./ -o ./ \
--use-bases-mask Y*,I*,I*,Y* \
--minimum-trimmed-read-length 35 \
--mask-short-adapter-reads 22 \
--sample-sheet ./SampleSheet_230915.csv \
--stats-dir ./Stats/ \
--reports-dir ./Reports/ \
--ignore-missing-bcls \

echo "Finish"
