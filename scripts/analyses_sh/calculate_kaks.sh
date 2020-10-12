#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20GB
#SBATCH --job-name calculate_KaKs
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
sample=$(pwd | sed s/^.*\\///)
SCRIPTS=$(pwd | sed s/angiosperms.*/angiosperms\\/scripts/)
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export PATH="$SCRIPTS/calculate_Ka_Ks_pipeline:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"
export PERL5LIB="$HOME/miniconda3/envs/gene-duplication/lib/perl5/site_perl/5.22.0"

#Calculate KaKs
cd dupgen/results
mkdir kaks_results
for i in *pairs
do
	echo "Calculating KaKs for each gene pair"
	pairs=$(echo $i | tr '.' '\t' | cut -f2)
	perl $SCRIPTS/calculate_Ka_Ks_pipeline/calculate_Ka_Ks_pipe.pl \
		-d ../../ref/mcscanx/"$sample"-cds.fa \
		-g "$i" \
		-o kaks_results/"$sample"."$pairs".kaks \
		-p $SCRIPTS/calculate_Ka_Ks_pipeline
	sed -i '2d' kaks_results/"$sample"."$pairs".kaks
done


