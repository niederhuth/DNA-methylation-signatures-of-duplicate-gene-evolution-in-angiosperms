#!/bin/bash --login
#SBATCH --time=1:00:00
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
outgroup=$(awk -v FS="," -v a="$sample" '$1 == a' ../../misc/genomes.csv | cut -d ',' -f 8)
collinear_genes=5
proximal_distance=10

#Copy over data
mkdir dupgen
mkdir dupgen/data
echo "Sample: $sample"
echo "Outgroup: $outgroup"
cp ref/mcscanx/"$sample".gff dupgen/data/"$sample".gff
cp diamond/"$sample"-"$sample".m8 dupgen/data/"$sample"_unfiltered.blast
cat ref/mcscanx/"$sample".gff ../"$outgroup"/ref/mcscanx/"$outgroup".gff > dupgen/data/"$sample"_"$outgroup".gff
cp diamond/"$sample"-"$outgroup".m8 dupgen/data/"$sample"_"$outgroup"_unfiltered.blast
cp ref/mcscanx/"$sample"_orthogroups.tsv dupgen/data/
cp ../"$outgroup"/ref/mcscanx/"$outgroup"_orthogroups.tsv dupgen/data/

#Filter BLAST hits based on orthogroups
cd dupgen/data
join -1 2 -2 1 <(sort -k2,2 "$sample"_unfiltered.blast) <(sort -k1,1 "$sample"_orthogroups.tsv) | \
tr ' ' '\t' > tmp
join -1 2 -2 1 <(sort -k2,2 tmp) <(sort -k1,1 "$sample"_orthogroups.tsv) | tr ' ' '\t' | \
awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | sort | uniq > "$sample".blast

join -1 2 -2 1 <(sort -k2,2 "$sample"_"$outgroup"_unfiltered.blast) <(sort -k1,1 "$outgroup"_orthogroups.tsv) | \
tr ' ' '\t' > tmp
join -1 2 -2 1 <(sort -k2,2 tmp) <(sort -k1,1 "$sample"_orthogroups.tsv) | tr ' ' '\t' | \
awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | sort | uniq > "$sample"_"$outgroup".blast
rm tmp

#Run DupGen_finder
cd ../
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
fgrep -w -v -f tmp "$sample".gff.sorted > tmp2
awk -v OFS="\t" '{ print $1,"unclassified" }' tmp2 >> classified_genes.tsv
rm tmp tmp2

cut -f2 classified_genes.tsv | sort | uniq -c

echo "Done"

