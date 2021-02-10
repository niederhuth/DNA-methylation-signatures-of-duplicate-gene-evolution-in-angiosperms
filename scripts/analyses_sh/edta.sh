#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80GB
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
fasta=../../"$sample"-edta.fa
cds="$sample"-cds.fa
gff="$sample"-edta.gff
bed="$sample"-cds.bed
threads=10
if [ $sample == "Zmays" ] 
then 
	species="Maize"
elif [ $sample == "Osativa" ]
then 
	species="Rice"
else 
	species="others"
fi

#Set up working directory
cd ref/annotations
mkdir edta
cp $cds edta/
#Create CDS bed file
echo "Creating CDS bed file"
awk '$3=="CDS"' $gff | convert2bed --input=gff - > edta/"$bed"
cd edta

#Run EDTA
echo "Running EDTA for $sample"
perl $HOME/miniconda3/envs/edta/bin/EDTA.pl \
	--genome $fasta \
	--species $species \
	--step all \
	--overwrite 1 \
	--cds $cds \
	--sensitive 0 \
	--anno 1 \
	--evaluate 0 \
	--exclude $bed \
	--force 1 \
	--threads $threads

