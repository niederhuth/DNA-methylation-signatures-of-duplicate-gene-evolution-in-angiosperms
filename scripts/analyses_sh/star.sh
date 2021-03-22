#!/bin/bash --login
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
#SBATCH --job-name align_rna
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda, do not change these
export PATH="${conda}/envs/gene-duplication/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/gene-duplication/lib:${LD_LIBRARY_PATH}"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*//)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
index="$(pwd | sed s/data.*/data/)/${species}/ref/STAR"
fastq="*.fastq.gz"

#Run Star
echo "Running STAR for ${sample}"
STAR \
	--runThreadN ${threads} \
	--runMode alignReads \
	--genomeDir ${index} \
	--readFilesIn ${fastq} \
	--outFileNamePrefix ${sample} \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMstrandField intronMotif \
	--outFilterType BySJout \
	--outFilterMultimapNmax 10 \
	--alignSJoverhangMin 5 \
	--alignSJDBoverhangMin 3 \
	--alignIntronMin 20 \
	--alignIntronMax 0 \
	--outFilterScoreMinOverLread 0.33 \
	--outFilterMatchNminOverLread 0.33 \
	--outFilterMismatchNmax 10 \
	--outFilterMismatchNoverReadLmax 0.1 \
	--quantMode GeneCounts

echo "Done"


