#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name=rescue_nonTE_orthogroups
#SBATCH --output=%x-%j.SLURMout

#Set variables
number_of_species=5 #min number of species in an orthogroup
percent_cutoff=0.34 #max percentage of species in an orthogroup with putative TE hits

orthogroup_size=5 #min number genes in orthogroup to consider for false positive calls
perc_genes=0.25 #max percentage of genes called in an orthogroup for it to be considered a false positive, set to 0 to turn off

#Count up orthogroups with putative TE-like genes and percentage of that orthogroup across species identified as TE-like
mkdir TE_filtering/
cat */ref/TE_filtering/orthogroups.txt | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > TE_filtering/tmp

for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		cut -f2 ${i}/ref/mcscanx/${i}_orthogroups.tsv | sort | uniq | awk -v OFS="\t" -v a=${i} '{print a,$0}' >> TE_filtering/tmp2
	fi
done
cd TE_filtering/
cut -f2 tmp2 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > tmp3
join -1 2 -2 2 tmp3 tmp | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$3/$2}' > filtered_ortholog_counts.tsv
rm tmp*

#Get a list of orthogroups to keep
awk -v a=${number_of_species} -v b=${percent_cutoff} '$2>=a && $4<=b' filtered_ortholog_counts.tsv | cut -f1 > orthogroups_to_keep.tsv

#Get a list of putative TE-like orthogroups
fgrep -v -f orthogroups_to_keep.tsv filtered_ortholog_counts.tsv | cut -f1 > putative_TE_orthogroups.tsv

#Get genes to keep for that species and annotations info if available
cd ../
for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		echo ${i}
		fgrep -f TE_filtering/orthogroups_to_keep.tsv ${i}/ref/TE_filtering/filtered_genes.tsv | sort -k1,1 > TE_filtering/${i}_genes_to_keep.tsv
		cut -f1 TE_filtering/${i}_genes_to_keep.tsv > TE_filtering/tmp
		if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
		then
			fgrep -f TE_filtering/tmp ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort -k1,1 > TE_filtering/tmp2
			join -1 1 -2 1 TE_filtering/${i}_genes_to_keep.tsv TE_filtering/tmp2 | tr ' ' '\t' | tr ';' ' ' > TE_filtering/tmp3
			mv TE_filtering/tmp3 TE_filtering/${i}_genes_to_keep.tsv
		fi
		cat ${i}/ref/TE_filtering/noTE_gene_list.txt TE_filtering/tmp > ${i}/ref/TE_filtering/noTE_gene_list_2.txt
		fgrep -v -f ${i}/ref/TE_filtering/noTE_gene_list_2.txt ${i}/ref/TE_filtering/filtered_genes.tsv > ${i}/ref/TE_filtering/filtered_genes_2.tsv
		rm TE_filtering/tmp*
	fi
done

#Get putative TE-like genes for that species and annotations info if available
for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		echo ${i}
		fgrep -f TE_filtering/putative_TE_orthogroups.tsv ${i}/ref/TE_filtering/filtered_genes.tsv | sort -k1,1 > TE_filtering/${i}_putative_TEs.tsv
		cut -f1 TE_filtering/${i}_putative_TEs.tsv > TE_filtering/tmp
		if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
		then
			fgrep -f TE_filtering/tmp ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort -k1,1 > TE_filtering/tmp2
			join -1 1 -2 1 TE_filtering/${i}_putative_TEs.tsv TE_filtering/tmp2 | tr ' ' '\t' | tr ';' ' ' > TE_filtering/tmp3
			mv TE_filtering/tmp3 TE_filtering/${i}_putative_TEs.tsv
		fi
		rm TE_filtering/tmp*
	fi
done





awk -v a=${orthogroup_size} -v b=${perc_genes} '$2>=a && $4<=b' filtered_orthgroup_counts.tsv | cut -f1 > good_orthogroups.txt
fgrep -f good_orthogroups.txt filtered_genes.tsv > genes_rescued_by_orthogroup.tsv
fgrep -v -f good_orthogroups.txt filtered_genes.tsv > filtered_genes_wo_rescued_TEs.tsv
cut -f1 genes_restored_by_orthogroup.tsv > tmp4
cat noTE_gene_list.txt tmp4 | sort > noTE_gene_list_w_rescued_genes.tsv
rm tmp*
