
for i in *
do
	if [ -d ${i}/ref/gene_filtering ]
	then
		echo $i
		cd ${i}/ref/gene_filtering
		cat pfam_filtered_genes.txt blast_filtered_genes.txt gypsy_filtered_genes.txt | sort | uniq > filtered_genes.tsv
		fgrep -f filtered_genes.tsv ../mcscanx/${i}_orthogroups.tsv | sort -k1,1 > tmp
		if [ -f ../../methylpy/results/${i}_classified_genes.tsv ]
		then
			fgrep -f filtered_genes.tsv ../../methylpy/results/${i}_classified_genes.tsv | cut -f1,30 | sort -k1,1 > tmp2
			cut -f1 tmp2 > tmp3
			fgrep -v -f tmp3 filtered_genes.tsv > tmp4
			awk -v OFS="\t" '{print $1,"NA"}' tmp4 >> tmp2
			sort -k1,1 tmp2 > tmp5
			join -1 1 -2 1 tmp tmp5 | tr ' ' '\t' > filtered_genes.tsv
		else
			mv tmp filtered_genes.tsv
		fi
		cut -f2 filtered_genes.tsv | sort | uniq > orthogroups.txt
		cut -f2 filtered_genes.tsv | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' | sort -k2,2 > tmp
		cut -f2 tmp > tmp2
		fgrep -f tmp2 ../mcscanx/${i}_orthogroups.tsv | cut -f2 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' | sort -k2,2 > tmp3
		join -1 2 -2 2 tmp3 tmp | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$3/$2}' > filtered_ortholog_counts.tsv
		rm tmp*
		cd ../../../
	fi
done



#All genes
mkdir orthofinder/gene_filtering/
cat */ref/gene_filtering/orthogroups.txt | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > orthofinder/gene_filtering/tmp

for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		cut -f2 ${i}/ref/mcscanx/${i}_orthogroups.tsv | sort | uniq | awk -v OFS="\t" -v a=${i} '{print a,$0}' >> orthofinder/gene_filtering/tmp2
	fi
done
cd orthofinder/gene_filtering/
cut -f2 tmp2 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > tmp3
join -1 2 -2 2 tmp3 tmp | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$3/$2}' > filtered_ortholog_counts.tsv
rm tmp*


awk '$2>40 && $4<0.5' filtered_ortholog_counts.tsv | cut -f1 > tmp5

fgrep -f tmp5 Athaliana/ref/gene_filtering/filtered_genes.tsv | cut -f3 | sort | uniq -c
