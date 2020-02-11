#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=48GB
#SBATCH --job-name recall_mC
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#List Variables
sample=$(pwd | sed s/^.*\\///)
fasta="../ref/$sample.fa"
unmethylated_control=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 4)

#Run Methylpy
echo "Unmethylated Control is $unmethylated_control"
cd methylpy
rm allc_"$sample".tsv.gz*

if ls fastq/*_2.fastq.gz >/dev/null 2>&1
then
	echo "Data is paired-end"
	paired_end="True"
else
	echo "Data is single-end"
	paired_end="False"
fi

echo "Recalling methylation states"
methylpy call-methylation-state \
	--input-file "$sample"_processed_reads_no_clonal.bam \
	--sample $sample \
	--ref-fasta $fasta \
	--paired-end $paired_end \
	--path-to-output "" \
	--num-procs 20 \
	--num-upstream-bases 0 \
	--num-downstream-bases 2 \
	--generate-allc-file True \
	--generate-mpileup-file True \
	--compress-output True \
	--bgzip False \
	--path-to-bgzip "" \
	--path-to-tabix "" \
	--path-to-samtools "" \
	--remove-chr-prefix False \
	--add-snp-info False \
	--unmethylated-control $unmethylated_control \
	--binom-test True \
	--sig-cutoff .01 \
	--min-mapq 30 \
	--min-cov 3 \
	--min-base-quality 1 \
	--keep-temp-files True

echo "Done"
