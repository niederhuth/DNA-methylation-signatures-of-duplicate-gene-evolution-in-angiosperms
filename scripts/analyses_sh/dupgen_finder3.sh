#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name DupGen_finder
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
SCRIPTS=$(pwd | sed s/angiosperms.*/angiosperms\\/scripts/)
export PATH="$HOME/miniconda3/envs/gene-duplication/bin:$PATH"
export PATH="$SCRIPTS/DupGen_finder:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#export PATH="../../../scripts/DupGen_finder:$PATH"

#Set Variables
sample=$(pwd | sed s/^.*\\///)
outgroup=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 6 | cut -d ' ' -f1)
collinear_genes=4
proximal_distance=10

#Copy over data
mkdir dupgen3
mkdir dupgen3/data
echo "Sample: $sample"
echo "Outgroup: $outgroup"
cp ref/mcscanx/"$sample".gff dupgen3/data/"$sample".gff
cp diamond/"$sample"-"$sample".m8 dupgen3/data/"$sample".blast
cat ref/mcscanx/"$sample".gff ../"$outgroup"/ref/mcscanx/"$outgroup".gff > dupgen3/data/"$sample"_"$outgroup".gff
cp diamond/"$sample"-"$outgroup".m8 dupgen3/data/"$sample"_"$outgroup".blast
cp ref/mcscanx/"$sample"_orthogroups.tsv dupgen3/data/
cp ../"$outgroup"/ref/mcscanx/"$outgroup"_orthogroups.tsv dupgen3/data/


#Run DupGen_finder
cd dupgen3
echo "Running DupGen_finder"
DupGen_finder.pl \
	-i data \
	-t $sample \
	-c $outgroup \
	-s $collinear_genes \
	-d $proximal_distance \
	-o results

echo "Running DupGen_finder-unique"
DupGen_finder-unique.pl \
        -i data \
        -t $sample \
        -c $outgroup \
	-s $collinear_genes \
	-d $proximal_distance \
        -o results-unique

echo "Make classified gene list"
cd results-unique
echo -e "Feature"'\t'"Duplication" > classified_genes.tsv
for i in *genes-unique *singletons
do
        a=$(echo $i | sed s/\.genes-unique// | sed s/.*\\.//)
        sed '1d' $i | awk -v b=$a -v OFS="\t" '{ print $1,b }' >> classified_genes.tsv
done
cut -f1 classified_genes.tsv > tmp
fgrep -v -f tmp "$sample".gff.sorted > tmp2
awk -v OFS="\t" '{ print $1,"unclassified" }' tmp2 >> classified_genes.tsv
rm tmp tmp2

#gene orientation
awk -v OFS="\t"  '{if ($3=="gene") print $9,$7}' ../../ref/annotations/"$sample".gff | \ 
sed s/^.*\=// > orientation.tsv

echo "Done"

