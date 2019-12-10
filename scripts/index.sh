#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --job-name index
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set tmp directories
export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

#List Variables
sample=$(pwd | sed s/^.*\\///)
methylC=$(grep $sample ../../misc/genomes.csv | cut -d ',' -f 5)
aligner="bowtie2"

#Setup Methylpy
if [ $methylC == "yes" ]
then
	echo "Building Methylpy index"
	cd ref
	methylpy build-reference \
		--input-files $sample.fa \
		--output-prefix $sample \
		--num-procs 10 \
		--aligner $aligner
else
	echo "No methylation data"
fi

