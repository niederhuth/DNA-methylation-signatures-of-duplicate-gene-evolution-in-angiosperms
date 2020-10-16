dispersed genes with collinear orthologs

PAV -> compare to outgroup & use to discuss fractionation




sed '1d' Duplications.tsv | while read line
do
a=$(echo $line | cut -f2)
b=$(echo $line | cut -f3)
echo $line | cut -f6,7 | tr ' ' '\n' | tr '\t' '\n' | gsed s/$/,/ > tmp
gsed s/^/$a,$b,/ tmp | tr ',' '\t' >> reformat.tsv
done

