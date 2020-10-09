#!/bin/bash --login
#SBATCH --time=98:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=160GB
#SBATCH --job-name methylpy
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#List Variables
sample=$(pwd | sed s/^.*\\///)
f_ref="../ref/'$sample'_f"
r_ref="../ref/'$sample'_r"
fasta="../ref/$sample.fa"
read_single="fastq/*fastq.gz"
read1="fastq/*_1.fastq.gz"
read2="fastq/*_2.fastq.gz"
unmethylated_control=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 4)
adaptor_single="AGATCGGAAGAGCACACGTCTG"
adaptor1="AGATCGGAAGAGCACACGTCTGAAC"
adaptor2="AGATCGGAAGAGCGTCGTGTAGGGA"
aligner="bowtie2"
aligner_options="--very-sensitive -X 1000"
picard="$HOME/miniconda3/envs/gene-duplication/share/picard-2.21.2-0"

#Run Methylpy
echo "Unmethylated Control is $unmethylated_control"
cd methylpy
if ls fastq/*_2.fastq.gz >/dev/null 2>&1
then
	echo "Data is paired-end"
	echo "Running methylpy"

	methylpy paired-end-pipeline \
	--read1-files $read1 \
	--read2-files $read2 \
	--sample $sample \
	--forward-ref $f_ref \
	--reverse-ref $r_ref \
	--ref-fasta $fasta \
	--libraries "libA" \
	--path-to-output "" \
	--pbat False \
	--check-dependency False \
	--num-procs 20 \
	--sort-mem 5G \
	--num-upstream-bases 0 \
	--num-downstream-bases 2 \
	--generate-allc-file True \
	--generate-mpileup-file True \
	--compress-output True \
	--bgzip False \
	--path-to-bgzip "" \
	--path-to-tabix "" \
	--trim-reads True \
	--path-to-cutadapt "" \
	--path-to-aligner "" \
	--aligner $aligner \
	--aligner-options "$aligner_options" \
	--merge-by-max-mapq True \
	--remove-clonal True \
	--path-to-picard $picard \
	--keep-clonal-stats True \
	--java-options "" \
	--path-to-samtools "" \
	--adapter-seq-read1 $adaptor1 \
	--adapter-seq-read2 $adaptor2 \
	--remove-chr-prefix False \
	--add-snp-info False \
	--unmethylated-control $unmethylated_control \
	--binom-test True \
	--sig-cutoff .01 \
	--min-mapq 30 \
	--min-cov 3 \
	--max-adapter-removal 1 \
	--overlap-length 3 \
	--error-rate 0.1 \
	--min-qual-score 10 \
	--min-read-len 30 \
	--min-base-quality 1 \
	--keep-temp-files True

else
	echo "Data is single-end"
	echo "Running Methylpy"

	methylpy single-end-pipeline \
        --read-files $read_single \
        --sample $sample \
        --forward-ref $f_ref \
        --reverse-ref $r_ref \
        --ref-fasta $fasta \
        --libraries "libA" \
        --path-to-output "" \
        --pbat False \
        --check-dependency False \
        --num-procs 20 \
        --sort-mem 5G \
        --num-upstream-bases 0 \
        --num-downstream-bases 2 \
        --generate-allc-file True \
        --generate-mpileup-file True \
        --compress-output True \
        --bgzip False \
        --path-to-bgzip "" \
        --path-to-tabix "" \
        --trim-reads True \
        --path-to-cutadapt "" \
        --path-to-aligner "" \
        --aligner $aligner \
        --aligner-options "$aligner_options" \
        --merge-by-max-mapq True \
        --remove-clonal True \
        --path-to-picard $picard \
        --keep-clonal-stats True \
        --java-options "" \
        --path-to-samtools "" \
        --adapter-seq $adaptor_single \
        --remove-chr-prefix False \
        --add-snp-info False \
        --unmethylated-control $unmethylated_control \
        --binom-test True \
        --sig-cutoff .01 \
        --min-mapq 30 \
        --min-cov 3 \
        --max-adapter-removal 1 \
        --overlap-length 3 \
        --error-rate 0.1 \
        --min-qual-score 10 \
        --min-read-len 30 \
        --min-base-quality 1 \
        --keep-temp-files True
fi

echo "Done"
