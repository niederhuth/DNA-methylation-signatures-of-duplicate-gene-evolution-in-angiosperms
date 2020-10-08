#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name DupGen_finder2
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
collinear_genes=4
proximal_distance=10

#Copy over data
mkdir dupgen2
mkdir dupgen2/data
echo "Sample: $sample"
echo "Outgroup: $outgroup"
cp ref/mcscanx/"$sample".gff dupgen2/data/"$sample".gff
cp diamond2/"$sample"-"$sample".m8 dupgen2/data/"$sample"_unfiltered.blast
cat ref/mcscanx/"$sample".gff ../"$outgroup"/ref/mcscanx/"$outgroup".gff > dupgen2/data/"$sample"_"$outgroup".gff
cp diamond2/"$sample"-"$outgroup".m8 dupgen2/data/"$sample"_"$outgroup"_unfiltered.blast
cp ref/mcscanx/"$sample"_orthogroups.tsv dupgen2/data/
cp ../"$outgroup"/ref/mcscanx/"$outgroup"_orthogroups.tsv dupgen2/data/

#Filter BLAST hits based on orthogroups
cd dupgen2/data
cat "$sample"_unfiltered.blast | while read line
do
        a=$(echo $line | cut -d ' ' -f1)
        b=$(echo $line | cut -d ' ' -f2)
        c=$(grep $a "$sample"_orthogroups.tsv | cut -f2)
        d=$(grep $b "$sample"_orthogroups.tsv | cut -f2)
        echo $line $c $d | tr ' ' '\t' | awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' >> "$sample".blast
done
cat "$sample"_"$outgroup"_unfiltered.blast | while read line
do
        a=$(echo $line | cut -d ' ' -f1)
        b=$(echo $line | cut -d ' ' -f2)
        c=$(grep $a "$sample"_orthogroups.tsv | cut -f2)
        d=$(grep $b "$outgroup"_orthogroups.tsv | cut -f2)
        echo $line $c $d | tr ' ' '\t' | awk -v OFS="\t" '{if ($13==$14) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' >> "$sample"_"$outgroup".blast
done

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
fgrep -v -f tmp "$sample".gff.sorted > tmp2
awk -v OFS="\t" '{ print $1,"unclassified" }' tmp2 >> classified_genes.tsv
rm tmp tmp2

#gene orientation
awk -v OFS="\t"  '{if ($3=="gene") print $9,$7}' ../../ref/annotations/"$sample".gff | \ 
sed s/^.*\=// > orientation.tsv

echo "Done"

