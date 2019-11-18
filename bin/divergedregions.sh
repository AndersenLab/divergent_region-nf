#!/bin/bash
filename='diverged_pops_red.tsv'
filelines=`cat $filename`
for line in $filelines ; do
	for chrom in I II III IV V X; do
	bcftools query -i 'GT="alt"' -r $chrom -s $line -f '[%CHROM\t%POS\t%GT\n]' Ce330_STRELKA_CSQ_annotated.vcf.gz |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E5  AS INT)*1E5, count(*) from T group by c,CAST(p/1E5  AS INT)*1E5 ; drop table T;");}' |\
	sqlite3 -separator ' ' $chrom.sqlite > $chrom.$line.txt
	rm $chrom.sqlite
    done;
done;