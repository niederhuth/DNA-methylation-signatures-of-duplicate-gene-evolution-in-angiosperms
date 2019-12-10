#!/bin/bash --login
#SBATCH --time=3:55:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name DupGen_finder
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
SCRIPTS=$(pwd | sed s/methylation.*/methylation\\/scripts/)
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export PATH="$SCRIPTS/DupGen_finder:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#export PATH="../../../scripts/DupGen_finder:$PATH"

#Set Variables
sample=$(pwd | sed s/^.*\\///)
outgroup=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 8)

#Run DupGen_finder
mkdir dupgen
mkdir dupgen/data
echo "Sample: $sample"
echo "Outgroup: $outgroup"
cp ref/mcscanx/"$sample".gff dupgen/data/"$sample".gff
cp diamond/"$sample"-"$sample".m8 dupgen/data/"$sample".blast
cat ref/mcscanx/"$sample".gff ../"$outgroup"/ref/mcscanx/"$outgroup".gff > dupgen/data/"$sample"_"$outgroup".gff
cp diamond/"$sample"-"$outgroup".m8 dupgen/data/"$sample"_"$outgroup".blast

cd dupgen
echo "Running DupGen_finder"
DupGen_finder.pl \
	-i data \
	-t $sample \
	-c $outgroup \
	-o results

echo "Running DupGen_finder-unique"
DupGen_finder-unique.pl \
        -i data \
        -t $sample \
        -c $outgroup \
        -o results-unique

echo "Make classified gene list"
cd results-unique
echo -e "Feature"'\t'"Duplication" >> classified_genes.tsv
for i in *genes-unique *singletons
do
        a=$(echo $i | sed s/\.genes-unique// | sed s/.*\\.//)
        sed '1d' $i | awk -v b=$a -v OFS="\t" '{ print $1,b }' >> classified_genes.tsv
done
cut -f1 classified_genes.tsv > tmp
fgrep -v -f tmp "$sample".gff.sorted > tmp2
awk -v OFS="\t" '{ print $1,"unclassified" }' tmp2 >> classified_genes.tsv
rm tmp tmp2

echo "Done"

