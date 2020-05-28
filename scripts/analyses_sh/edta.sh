#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=60GB
#SBATCH --job-name edta
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/edta/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/edta/lib:$LD_LIBRARY_PATH"

#Set tmp directories
export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

#variables
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)
fasta=../"$sample".fa
cds="$sample"-cds.fa
gff="$sample"-cds.gff
bed="$sample"-cds.bed
threads=
if [ $sample == "Zmays" ] 
then 
	species="Maize"
elif [ $sample == "Osativa" ]
then 
	species="Rice"
else 
	species="others"
fi

#Create CDS bed file
cd ref/annotations
echo "Creating CDS bed file"
awk '$3=="CDS"' $gff | convert2bed --input=gff - > $bed
#Run EDTA
echo "Running EDTA for $sample"
perl $HOME/miniconda3/envs/edta/bin/EDTA.pl
	--genome $fasta
	--species $species
	--step all
	--overwrite	1
	--cds $cds
	--sensitive 0
	--anno 1
	--evaluate 0
	--exclude $bed
	--threads $threads

