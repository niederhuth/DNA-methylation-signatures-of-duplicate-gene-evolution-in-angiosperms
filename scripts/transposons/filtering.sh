#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name=good_orthogroups
#SBATCH --output=%x-%j.SLURMout

#Set variables
number_of_species=5
percent_cutoff=0.5

#Cound up orthogroups with putative TE-like genes and percentage of that orthogroup across species identified as TE-like
mkdir gene_filtering/
cat */ref/gene_filtering/orthogroups.txt | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > gene_filtering/tmp

for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		cut -f2 ${i}/ref/mcscanx/${i}_orthogroups.tsv | sort | uniq | awk -v OFS="\t" -v a=${i} '{print a,$0}' >> gene_filtering/tmp2
	fi
done
cd gene_filtering/
cut -f2 tmp2 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > tmp3
join -1 2 -2 2 tmp3 tmp | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$3/$2}' > filtered_ortholog_counts.tsv
rm tmp*

#Get a list of orthogroups to keep
awk -v a=${number_of_species} -v b=${percent_cutoff} '$2>a && $4<b' filtered_ortholog_counts.tsv | cut -f1 > orthogroups_to_keep.tsv

#Get a list of putative TE-like orthogroups
fgrep -v -f orthogroups_to_keep.tsv filtered_ortholog_counts.tsv | cut -f1 > putative_TE_orthogroups.tsv

#Get genes to keep for that species and annotations info if available
cd ../
for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		echo ${i}
		fgrep -f gene_filtering/orthogroups_to_keep.tsv ${i}/ref/gene_filtering/filtered_genes.tsv | sort -k1,1 > gene_filtering/${i}_genes_to_keep.tsv
		cut -f1 gene_filtering/${i}_genes_to_keep.tsv > gene_filtering/tmp
		if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
		then
			fgrep -f gene_filtering/tmp ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort -k1,1 > gene_filtering/tmp2
			join -1 1 -2 1 gene_filtering/${i}_genes_to_keep.tsv gene_filtering/tmp2 | tr ' ' '\t' | tr ';' ' ' > gene_filtering/tmp3
			mv gene_filtering/tmp3 gene_filtering/${i}_genes_to_keep.tsv
		fi
		cat ${i}/ref/gene_filtering/noTE_gene_list.txt gene_filtering/tmp > ${i}/ref/gene_filtering/noTE_gene_list_2.txt
		fgrep -v -f ${i}/ref/gene_filtering/noTE_gene_list_2.txt ${i}/ref/gene_filtering/filtered_genes.tsv > ${i}/ref/gene_filtering/filtered_genes_2.tsv
		rm gene_filtering/tmp*
	fi
done

#Get putative TE-like genes for that species and annotations info if available
for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		echo ${i}
		fgrep -f gene_filtering/putative_TE_orthogroups.tsv ${i}/ref/gene_filtering/filtered_genes.tsv | sort -k1,1 > gene_filtering/${i}_putative_TEs.tsv
		cut -f1 gene_filtering/${i}_putative_TEs.tsv > gene_filtering/tmp
		if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
		then
			fgrep -f gene_filtering/tmp ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort -k1,1 > gene_filtering/tmp2
			join -1 1 -2 1 gene_filtering/${i}_putative_TEs.tsv gene_filtering/tmp2 | tr ' ' '\t' | tr ';' ' ' > gene_filtering/tmp3
			mv gene_filtering/tmp3 gene_filtering/${i}_putative_TEs.tsv
		fi
		rm gene_filtering/tmp*
	fi
done
