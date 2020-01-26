#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --job-name setup
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
methylC=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 5)
aligner="bowtie2"

#Build Diamond database
cd ref/mcscanx
echo "Building Diamond database"
diamond makedb \
	--in "$sample"-protein.fa \
	--db "$sample".dmnd
cd ../../

#Setup Methylpy
if [ $methylC == "yes" ]
then
	mkdir methylpy
	echo "Building Methylpy index"
	cd ref
	methylpy build-reference \
		--input-files $sample.fa \
		--output-prefix $sample \
		--num-procs 4 \
		--aligner $aligner
	cd ../methylpy
	mkdir fastq results
	cd fastq
	echo "Downloading SRA files"
	sra_list=$(grep "$sample" ../../../../misc/sra_list.txt | cut -f 2 | tr ',' ' ')
        for i in $sra_list
        do
		echo "Downloading $i"
		prefetch --max-size 100000000 $i
		mv $i/*sra ./
		rmdir $i
                echo "Converting $i"
                fastq-dump --gzip --split-3 "$i".sra
		#gzip "$i"*fastq
                rm "$i".sra
        done
else
	echo "No methylation data"
fi

