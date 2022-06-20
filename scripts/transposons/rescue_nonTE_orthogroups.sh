#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name=rescue_nonTE_orthogroups
#SBATCH --output=%x-%j.SLURMout

#Set variables
perc_species_filter=TRUE #Filter based on percentage of species in an orthogroup with TE-like hits; TRUE or FALSE
number_of_species=5 #Min number of species in an orthogroup
species_perc_cutoff=0.2 #Max percentage of species in an orthogroup with putative TE hits
phylo_gene_filter=TRUE #Filter based on percentage of genes in across phylogeny with TE-like hits; TRUE or FALSE
phylo_number_species=2 #Min number of species an orthogroups has to be present in to be considered
phylo_number_genes=5 #Min number of genes in an orthogroup across the phylogeny to consider the orthogroup
phylo_perc_gene_cutoff=0.2 #Max percentage of genes in an orthogroup across the phylogeny to retain the orthogroup
species_gene_filter=TRUE #Filter based on percentage of genes within a species with TE-like hits; TRUE or FALSe
species_number_genes=5 #Min number of genes in an orthogroup within the species to consider the orthogroup
species_perc_gene_cutoff=0.2 #Max percentage of genes in an orthogroup within the species to retain the orthogroup

#Change to current directory
cd ${PBS_O_WORKDIR}
#The following shouldn't need to be changed, but should set automatically
path1="TE_filtering"

#Check for output directory, if not, make it
if [ -d ${path1} ]
then 
	echo "${path1} already exists, if you wish to redo this analysis, delete ${path1} and resubmit"
	exit 0
else
	mkdir ${path1}
fi

#Count up orthogroups with putative TE-like genes and percentage of that orthogroup across species identified as TE-like
#Combine list of filtered orthogroups and their individual counts
cat */ref/${path1}/filtered_orthgroup_counts.tsv > ${path1}/tmp
#Combine the uniq orthogroups for each species and output a table with two columns, column 1 is speces, column 2 the orthogroup
for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		cut -f2 ${i}/ref/mcscanx/${i}_orthogroups.tsv | sort | uniq | awk -v OFS="\t" -v a=${i} '{print a,$0}' >> ${path1}/species_orthogroup_list.txt
	fi
done
cd ${path1}/
#We are going to count up all the genes and all the putative TE-like genes across all species for each filtered orthogroup
#Output is: column 1: orthogroup; column 2: total gene count across all species; column 3: TE-like gene count across all species;
#column 4: percentage of TE-like genes across all species
cut -f1 tmp | sort | uniq | while read line
do
	grep ${line} tmp | awk -v OFS="\t" '{a+=$2;b+=$3}END{print $1,a,b,b/a}' >> tmp2
done
sort -k1,1 tmp2 > filtered_orthogroup_gene_counts.tsv
#Count up number of species per orthogroup
cut -f2 species_orthogroup_list.txt | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > tmp3
#Count up the number of species per filtered orthogroup
cut -f1 tmp | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > tmp4
join -1 2 -2 2 tmp3 tmp4 | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$3/$2}' | sort -k1,1 > filtered_orthogroup_counts.tsv
rm tmp*
#Combine some files
#Add a header
echo "Orthogroup Total_Species Total_Species_TE Perc_Species_TE Total_Genes Total_Genes_TE Perc_Species_TE" | tr ' ' '\t' > tmp
join -1 1 -2 1 filtered_orthogroup_counts.tsv filtered_orthogroup_gene_counts.tsv | tr ' ' '\t'> tmp2
cat tmp tmp2 > filtered_orthogroups.tsv
rm tmp tmp2

#Apply various filters to get orthogroups with likely false-positive genes that we want to keep
#There are two possible criteria we can use to rescue orthogroups at the phylogeny level
#Filter based on number of species with a TE-like hit
if [ perc_species_filter="TRUE" ]
then
	echo "Filtering based on the percentage of species with filtered as TE-like"
	awk -v a=${number_of_species} -v b=${species_perc_cutoff} '$2>=a && $4<=b' filtered_orthogroups.tsv > perc_species_filter_keep.tsv
	cut -f1 perc_species_filter_keep.tsv > perc_species_filter_keep_list.txt
	mkdir perc_species
fi
#Filter based on percentage of genes across the phylogeny
if [ phylo_gene_filter="TRUE" ]
then
	echo "Filtering based on the percentage of genes across the phylogeny with filtered as TE-like"
	awk -v a=${phylo_number_species} -v b=${phylo_number_genes} \
	-v c=${phylo_perc_gene_cutoff} '$2>=a && $5>=b && $7<=c' filtered_orthogroups.tsv > phylo_genes_filter_keep.tsv
	cut -f1 perc_species_filter_keep.tsv > phylo_genes_filter_keep_list.txt
	mkdir phylo_genes
fi
#For species filter, we will only create the directory, as this has to be applied to each individual species
if [ species_gene_filter="TRUE" ]
then
	mkdir species_genes
fi

#Get a list of putative TE-like orthogroups
#fgrep -v -f orthogroups_to_keep.tsv filtered_orthogroup_counts.tsv | cut -f1 > putative_TE_orthogroups.tsv

#Now loop over each species and get the rescued genes
mkdir putative_TEs rescued_genes
cd ../
for i in *
do
	if [ -f ${i}/ref/mcscanx/${i}_orthogroups.tsv ]
	then
		echo ${i}
		#Get data for prec_species_filter
		if [ perc_species_filter="TRUE" ]
		then
			#Get the list of genes for that species, with methylation info
			fgrep -f ${path1}/perc_species_filter_keep_list.txt ${i}/ref/${path1}/filtered_genes.tsv | \
			sort -k1,1 > ${path1}/perc_species/${i}_genes_to_keep.tsv
			#Get gene list
			cut -f1 ${path1}/perc_species/${i}_genes_to_keep.tsv > ${path1}/perc_species/tmp
			#Add annotation description if available
			if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
			then
				#Grep the genes
				while read line
				do
					anno=$(grep ${line} ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort | uniq | cut -f2 | tr '\n' ';')
					echo "${line} ${anno}" | tr ' ' '\t' >> ${path1}/perc_species/tmp2
				done < ${path1}/perc_species/tmp
				sort -k1,1 ${path1}/perc_species/tmp2 > ${path1}/perc_species/tmp3
				#Join the tables
				join -1 1 -2 1 ${path1}/perc_species/${i}_genes_to_keep.tsv ${path1}/perc_species/tmp3 | \
				tr ' ' '\t' | tr ';' ' ' > ${path1}/perc_species/tmp4
				#Check if any genes were
				cut -f1 ${path1}/perc_species/tmp4 > ${path1}/perc_species/tmp5
				fgrep -v -f ${path1}/perc_species/tmp5 ${path1}/perc_species/${i}_genes_to_keep.tsv > ${path1}/perc_species/tmp6
				#Combine
				cat ${path1}/perc_species/tmp4 ${path1}/perc_species/tmp6 > ${path1}/perc_species/${i}_genes_to_keep.tsv
			fi
			#create a list of genes to rescue
			awk -v OFS="\t" '{print $0,"perc_species"}' ${path1}/perc_species/${i}_genes_to_keep.tsv >> ${path1}/rescued_genes/${i}_rescue_genes.tsv
			#cleanup
			rm ${path1}/perc_species/tmp*
		fi
		#Get data for phylo_genes_filter
		if [ phylo_genes_filter="TRUE" ]
		then
			#Get the list of genes for that species, with methylation info
			fgrep -f ${path1}/phylo_genes_filter_keep_list.txt ${i}/ref/${path1}/filtered_genes.tsv | \
			sort -k1,1 > ${path1}/phylo_genes/${i}_genes_to_keep.tsv
			#Get gene list
			cut -f1 ${path1}/phylo_genes/${i}_genes_to_keep.tsv > ${path1}/phylo_genes/tmp
			#Add annotation description if available
			if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
			then
				#Grep the genes
				cat ${path1}/phylo_genes/tmp | while read line
				do
					anno=$(grep ${line} ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort | uniq | cut -f2 | tr '\n' ';')
					echo "${line} ${anno}" | tr ' ' '\t' >> ${path1}/phylo_genes/tmp2
				done
				sort -k1,1 ${path1}/phylo_genes/tmp2 > ${path1}/phylo_genes/tmp3
				#Join the tables
				join -1 1 -2 1 ${path1}/phylo_genes/${i}_genes_to_keep.tsv ${path1}/phylo_genes/tmp3 | \
				tr ' ' '\t' | tr ';' ' ' > ${path1}/phylo_genes/tmp4
				#Check if any genes were lost
				cut -f1 ${path1}/phylo_genes/tmp4 > ${path1}/phylo_genes/tmp5
				fgrep -v -f ${path1}/phylo_genes/tmp5 ${path1}/phylo_genes/${i}_genes_to_keep.tsv > ${path1}/phylo_genes/tmp6
				#Combine
				cat ${path1}/phylo_genes/tmp4 ${path1}/phylo_genes/tmp6 > ${path1}/phylo_genes/${i}_genes_to_keep.tsv
			fi
			#create/add to a list of genes to rescue
			awk -v OFS="\t" '{print $0,"phylo_genes"}' ${path1}/phylo_genes/${i}_genes_to_keep.tsv >> ${path1}/rescued_genes/${i}_rescue_genes.tsv
			#cleanup
			rm ${path1}/phylo_genes/tmp*
		fi
		#Get data for species_genes_filter
		if [ species_genes_filter="TRUE" ]
		then
			#Apply the filter to the orthogroups for that species
			awk -v a=${species_number_genes} -v b=${species_perc_gene_cutoff} '$2>=a && $4<=b' $i/ref/${path1}/filtered_orthgroup_counts.tsv | 
			cut -f1 > ${path1}/species_genes/tmp
			#Get the list of genes for that species, with methylation info
			fgrep -f ${path1}/species_genes/tmp ${i}/ref/${path1}/filtered_genes.tsv | \
			sort -k1,1 > ${path1}/species_genes/${i}_genes_to_keep.tsv
			cut -f1 ${path1}/species_genes/${i}_genes_to_keep.tsv > ${path1}/species_genes/tmp2
			#Add annotation description if available
			if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
			then
				#Grep the genes
				cat ${path1}/species_genes/tmp2 | while read line
				do
					anno=$(grep ${line} ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort | uniq | cut -f2 | tr '\n' ';')
					echo "${line} ${anno}" | tr ' ' '\t' >> ${path1}/species_genes/tmp3
				done
				sort -k1,1 ${path1}/species_genes/tmp3 > ${path1}/species_genes/tmp4
				#Join the tables
				join -1 1 -2 1 ${path1}/species_genes/${i}_genes_to_keep.tsv ${path1}/species_genes/tmp4 | \
				tr ' ' '\t' | tr ';' ' ' > ${path1}/species_genes/tmp5
				#Check if any genes were lost
				cut -f1 ${path1}/species_genes/tmp5 > ${path1}/species_genes/tmp6
				fgrep -v -f ${path1}/species_genes/tmp6 ${path1}/species_genes/${i}_genes_to_keep.tsv > ${path1}/species_genes/tmp7
				#Combine
				cat ${path1}/species_genes/tmp5 ${path1}/species_genes/tmp7 > ${path1}/species_genes/${i}_genes_to_keep.tsv
			fi
			#create/add to a list of genes to rescue
			awk -v OFS="\t" '{print $0,"species_genes"}' ${path1}/species_genes/${i}_genes_to_keep.tsv >> ${path1}/rescued_genes/${i}_rescue_genes.tsv
			#cleanup
			rm ${path1}/species_genes/tmp*
		fi
		#Get temp list of rescued genes and add to noTE gene list
		cut -f1 ${path1}/rescued_genes/${i}_rescue_genes.tsv | sort | uniq > ${path1}/rescued_genes/tmp
		cat ${i}/ref/${path1}/noTE_gene_list.txt ${path1}/rescued_genes/tmp > ${i}/ref/${path1}/noTE_gene_list_w_rescued_genes.txt
		#Create list of putative TEs remaining
		fgrep -v -f ${path1}/rescued_genes/tmp ${i}/ref/${path1}/filtered_genes.tsv > ${path1}/putative_TEs/${i}_putative_TEs.tsv
		cut -f1 ${path1}/putative_TEs/${i}_putative_TEs.tsv > ${path1}/putative_TEs/tmp
		#Add annotation data if present for putative TEs
		if [ -f ${i}/ref/annotations/${i}-annotations.txt ]
			then
				#Grep the genes
				cat ${path1}/putative_TEs/tmp | while read line
				do
					anno=$(grep ${line} ${i}/ref/annotations/${i}-annotations.txt | cut -f2,13 | tr ' ' ';' | sort | uniq | cut -f2 | tr '\n' ';')
					echo "${line} ${anno}" | tr ' ' '\t' >> ${path1}/putative_TEs/tmp2
				done
				sort -k1,1 ${path1}/putative_TEs/tmp2 > ${path1}/putative_TEs/tmp3
				#Join the tables
				join -1 1 -2 1 ${path1}/putative_TEs/${i}_putative_TEs.tsv ${path1}/putative_TEs/tmp3 | \
				tr ' ' '\t' | tr ';' ' ' > ${path1}/putative_TEs/tmp4
				#Check if any genes were lost
				cut -f1 ${path1}/putative_TEs/tmp4 > ${path1}/putative_TEs/tmp5
				fgrep -v -f ${path1}/putative_TEs/tmp5 ${path1}/putative_TEs/${i}_putative_TEs.tsv > ${path1}/putative_TEs/tmp6
				#Combine
				cat ${path1}/putative_TEs/tmp4 ${path1}/putative_TEs/tmp6 > ${path1}/putative_TEs/${i}_putative_TEs.tsv
			fi
		#Copy list of putative TEs to species folder
		cp ${path1}/putative_TEs/${i}_putative_TEs.tsv ${i}/ref/${path1}/${i}_putative_TEs.tsv
		#Cleanup
		rm ${path1}/rescued_genes/tmp* ${path1}/putative_TEs/tmp*
	fi
done

echo "Done"
