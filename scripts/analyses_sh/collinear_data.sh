rm collinear.bed tmp*
grep -v \# Athaliana.collinearity | sed 's/^ *//' | sed 's/- */-/' | sed 's/\://' | tr '-' '\t' | cut -f1-4 > tmp

for i in $(cut -f1 tmp | sort | uniq)
do
	awk -v a=$i '$1==a' tmp > tmp2
	genes=$(wc -l tmp2 | cut -d ' ' -f1)
	fgrep -f <(cut -f3 tmp2) ../../methylpy/results/Athaliana_classified_genes.tsv | cut -f23 | sort | uniq -c | sed 's/^ *//' > tmp3
	gbM1=$(grep gbM tmp3 | cut -d ' ' -f1)
	if [ -z "$gbM1" ]
	then
		gbM1="NA"
	fi
	TE1=$(grep TE-like tmp3 | cut -d ' ' -f1)
	if [ -z "$TE1" ]
	then
		TE1="NA"
	fi
	UM1=$(grep Unmethylated tmp3 | cut -d ' ' -f1)
	if [ -z "$UM1" ]
	then
		UM1="NA"
	fi
	UC1=$(grep Unclassified tmp3 | cut -d ' ' -f1)
	if [ -z "$UC1" ]
	then
		UC1="NA"
	fi
	fgrep -f <(cut -f4 tmp2) ../../methylpy/results/Athaliana_classified_genes.tsv | cut -f23 | sort | uniq -c | sed 's/^ *//' > tmp3
	gbM2=$(grep gbM tmp3 | cut -d ' ' -f1)
	if [ -z "$gbM2" ]
	then
		gbM2="NA"
	fi
	TE2=$(grep TE-like tmp3 | cut -d ' ' -f1)
	if [ -z "$TE2" ]
	then
		TE2="NA"
	fi
	UM2=$(grep Unmethylated tmp3 | cut -d ' ' -f1)
	if [ -z "$UM2" ]
	then
		UM2="NA"
	fi
	UC2=$(grep Unclassified tmp3 | cut -d ' ' -f1)
	if [ -z "$UC2" ]
	then
		UC2="NA"
	fi
	chr1=$(grep $(head -1 tmp2 | cut -f3) ../../ref/mcscanx/Athaliana.gff | cut -f1 | sed s/Athaliana_//)
	chr2=$(grep $(head -1 tmp2 | cut -f4) ../../ref/mcscanx/Athaliana.gff | cut -f1 | sed s/Athaliana_//)
	start1=$(grep $(head -1 tmp2 | cut -f3) ../../ref/mcscanx/Athaliana.gff | cut -f3)
	start2=$(grep $(head -1 tmp2 | cut -f4) ../../ref/mcscanx/Athaliana.gff | cut -f3)
	stop1=$(grep $(tail -1 tmp2 | cut -f3) ../../ref/mcscanx/Athaliana.gff | cut -f4)
	stop2=$(grep $(tail -1 tmp2 | cut -f4) ../../ref/mcscanx/Athaliana.gff | cut -f4)
	if [ "$start1" -gt "$stop1" ]
	then
		echo $chr1 $stop1 $start1 "$i"-1 $genes $gbM1 $TE1 $UM1 $UC1 | tr ' ' '\t' >> collinear.bed
	else 
		echo $chr1 $start1 $stop1 "$i"-1 $genes $gbM1 $TE1 $UM1 $UC1 | tr ' ' '\t' >> collinear.bed
	fi
	if [ "$start2" -gt "$stop2" ]
	then
		echo $chr2 $stop2 $start2 "$i"-2 $genes $gbM2 $TE2 $UM2 $UC2 | tr ' ' '\t' >> collinear.bed
	else 
		echo $chr2 $start2 $stop2 "$i"-2 $genes $gbM2 $TE2 $UM2 $UC2 | tr ' ' '\t' >> collinear.bed
	fi
done

sed 's/^Athaliana_//' ../../ref/mcscanx/Athaliana.gff | awk -v OFS="\t" '{print $1,$3,$4,$2}' > tmp
bedtools intersect -a collinear.bed -b tmp > tmp2

rm collinear2.bed
for i in $(cut -f4 tmp2 | sort | uniq)
do
	awk -v a=$i '$4==a' tmp2 > tmp3
	line=$(awk -v a=$i '$4==a' collinear.bed)
	genes=$(wc -l tmp3 | cut -d ' ' -f1)
	fgrep -f <(cut -f13 tmp3) ../../methylpy/results/Athaliana_classified_genes.tsv | cut -f23 | sort | uniq -c | sed 's/^ *//' > tmp4
	gbM=$(grep gbM tmp4 | cut -d ' ' -f1)
	if [ -z "$gbM" ]
	then
		gbM="NA"
	fi
	TE=$(grep TE-like tmp4 | cut -d ' ' -f1)
	if [ -z "$TE" ]
	then
		TE="NA"
	fi
	UM=$(grep Unmethylated tmp4 | cut -d ' ' -f1)
	if [ -z "$UM" ]
	then
		UM="NA"
	fi
	UC=$(grep Unclassified tmp4 | cut -d ' ' -f1)
	if [ -z "$UC" ]
	then
		UC="NA"
	fi
	echo $line $genes $gbM $TE $UM $UC | tr ' ' '\t' >> collinear2.bed
done

echo "Chr Start Stop Name Collinear_Genes Collinear_gbM Collinear_TE.like Collinear_Unmethylated Collinear_Unclassified Total_Genes Total_gbM Total_TE.like Total_Unmethylated Total_Unclassified TE_Number TE_bps Total_Length Percent_TE" | tr ' ' '\t' > collinear_data.tsv
bedtools coverage -a collinear2.bed -b ../../ref/annotations/Athaliana-TEanno.gff >> collinear_data.tsv







