#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --job-name star-index
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=10
read_length=50

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/rna/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/rna/lib:${LD_LIBRARY_PATH}"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*// | sed s/$/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)

#Set the sjdbOverhang
#This should be the input read_length minus 1
sjdbOverhang=$(expr ${read_length} - 1)

#Change to ref directory
cd ref
#Get version
version=$(ls ${i}-v*.fa | sed s/.*\-v// | sed s/.fa//)

#Make gtf file from gff
echo "Converting gff to gtf"
gffread annotations/${genotype}-v${version}.gff -T > annotations/${genotype}-v${version}.gtf

#Create STAR index files
echo "Making STAR index"
STAR \
	--runThreadN ${threads} \
	--runMode genomeGenerate \
	--genomeDir STAR/ \
	--genomeFastaFiles ${species}.fa \
	--sjdbGTFfile annotations/${species}.gtf \
	--sjdbOverhang ${sjdbOverhang}

echo "Done"