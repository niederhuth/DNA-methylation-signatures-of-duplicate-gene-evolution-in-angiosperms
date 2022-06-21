#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name=../job_reports/filter_TEs
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/gene-duplication/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/gene-duplication/lib:$LD_LIBRARY_PATH"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
path3="TE_filtering"

#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#Generate a gene list from gff file
echo "Generating Gene List"
perl ${path2}/transposons/pl/create_gene_list.pl \
	--input_gff ../annotations/${species}.gff \
	--output_file input_gene_list.txt

#Check for Pfam A
echo "Downloading Pfam-A"
wget -q http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

#Prepare Pfam A hmm-database
echo "Preparing Pfam A hmm database"
hmmpress Pfam-A.hmm

#Search Pfam A hmm domains
echo "Searching against Pfam A hmm profiles"
hmmscan \
	--domE 1e-5 \
	-E 1e-5 \
	--cpu ${threads} \
	-o pfam_alignments.out \
	--tblout prot_domains.out \
	Pfam-A.hmm \
	../mcscanx/${species}-protein.fa

#Download and make Transposase blast DB
echo "Downloading Tpases020812"
wget -q http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz
gunzip Tpases020812.gz
#remove trailing white space in the Tpases020812 file
sed -i 's/[ \t]*$//' Tpases020812

echo "Making Transposase diamond blast DB"
diamond makedb \
	--threads ${threads} \
	--in Tpases020812 \
	--db Tpases020812.dmnd

#Run diamond blastp against Transposases 
echo "Running diamond blastp on Transposases"
diamond blastp \
	--threads ${threads} \
	--db Tpases020812.dmnd \
	--query ../mcscanx/${species}-protein.fa \
	--out TE_blast.out\
	--evalue 1e-10 \
	--outfmt 6

#Create a genelist with no TEs
echo "Creating gene list with TEs removed"
python ${path2}/transposons/py/create_no_TE_genelist.py \
	--input_geneList input_gene_list.txt \
	--pfamhmm prot_domains.out \
	--TEpfam_list ${path1}/TE_Pfam_domains.txt \
	--TEblast TE_blast.out \
	--output_file noTE_gene_list.txt

#Now combine with ortholog & methylation calls
cat pfam_filtered_genes.txt blast_filtered_genes.txt gypsy_filtered_genes.txt | sort | uniq > filtered_genes.tsv
fgrep -w -f filtered_genes.tsv ../mcscanx/${species}_orthogroups.tsv | sort -k1,1 > tmp
if [ -f ../../methylpy/results/${species}_classified_genes.tsv ]
then
	fgrep -w -f filtered_genes.tsv ../../methylpy/results/${species}_classified_genes.tsv | cut -f1,30 | sort -k1,1 > tmp2
	cut -f1 tmp2 > tmp3
	fgrep -w -v -f tmp3 filtered_genes.tsv > tmp4
	awk -v OFS="\t" '{print $1,"NA"}' tmp4 >> tmp2
	sort -k1,1 tmp2 > tmp5
	join -1 1 -2 1 tmp tmp5 | tr ' ' '\t' > filtered_genes.tsv
else
	mv tmp filtered_genes.tsv
fi

#Count the number & percentage of genes for filtered orthogroups
#Useful for reducing false positives
cut -f2 filtered_genes.tsv | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' | sort -k2,2 > tmp
cut -f2 tmp > tmp2
fgrep -w -f tmp2 ../mcscanx/${species}_orthogroups.tsv | cut -f2 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' | sort -k2,2 > tmp3
join -1 2 -2 2 tmp3 tmp | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$3/$2}' > filtered_orthgroup_counts.tsv
rm tmp*

echo "Done"

