#!/bin/bash
for wsz in 1E4 5E4 ; do
	for chrom in I II III IV V X; do
	bcftools query -i 'GT="alt"' -r $chrom -s CB4856 -f '[%CHROM\t%POS\t%GT\n]' Ce330_STRELKA_CSQ_annotated.vcf.gz |\
	awk -v ws="$wsz" -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E5  AS INT)*1E5, count(*) from T group by c,CAST(p/1E5  AS INT)*1E5 ; drop table T;");}' |\
	sqlite3 -separator ' ' $chrom.sqlite > $chrom.$line.txt
	rm $chrom.sqlite
    done;
done;

for wsz in 1E4 5E4 ; do
for chrom in I II III IV V X; do
	bcftools query -i 'GT="alt"' -r $chrom -s XZ1516 -f '[%CHROM\t%POS\t%GT\n]' Ce330_STRELKA_CSQ_annotated.vcf.gz |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CB.sqlite > $chrom.1E4_XZ1516.txt
	rm CB.sqlite
done;
done;

for chrom in I II III IV V X; do
	bcftools query -i 'GT="alt"' -r $chrom -s XZ1516 -f '[%CHROM\t%POS\t%GT\n]' Ce330_STRELKA_CSQ_annotated.vcf.gz |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CB.sqlite > $chrom.1E4_XZ1516.txt
	rm CB.sqlite
done;

for chrom in I II III IV V X; do
	bcftools query -i 'GT="alt"' -r $chrom -s ECA701 -f '[%CHROM\t%POS\t%GT\n]' Ce330_STRELKA_CSQ_annotated.vcf.gz |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CB.sqlite > $chrom.1E4_ECA701.txt
	rm CB.sqlite
done;

for chrom in I II III IV V X; do
	bcftools query -i 'GT="alt"' -r $chrom -s ECA732 -f '[%CHROM\t%POS\t%GT\n]' Ce330_STRELKA_CSQ_annotated.vcf.gz |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CB.sqlite > $chrom.1E4_ECA732.txt
	rm CB.sqlite
done;

for chrom in I II III IV V X; do
	bcftools query -i 'GT="alt"' -r $chrom -s DL238 -f '[%CHROM\t%POS\t%GT\n]' Ce330_STRELKA_CSQ_annotated.vcf.gz |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CB.sqlite > $chrom.1E4_DL238.txt
	rm CB.sqlite
done;


for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../input_files/Ce330_GATK4_STRELKA2_Intersect.vcf.gz |\
	bcftools query -i 'GT="alt"' -s XZ1516 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_XZ1516_GATK4.txt
	rm CBbe.sqlite
done;


for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../input_files/Ce330_GATK4_STRELKA2_Intersect.vcf.gz |\
	bcftools query -i 'GT="alt"' -s ECA701 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_ECA701_GATK4.txt
	rm CBbe.sqlite
done;

for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../input_files/Ce330_GATK4_STRELKA2_Intersect.vcf.gz |\
	bcftools query -i 'GT="alt"' -s CB4856 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_CB4856_GATK4.txt
	rm CBbe.sqlite
done;

for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../input_files/Ce330_GATK4_STRELKA2_Intersect.vcf.gz |\
	bcftools query -i 'GT="alt"' -s N2 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_N2_GATK4.txt
	rm CBbe.sqlite
done;

for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../20181113-GATK4/ANNOTATE_VCF/Ce330_annotated.vcf.gz |\
	bcftools query -i 'GT="alt"' -s N2 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_N2_GATK4_masked.txt
	rm CBbe.sqlite
done;

for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../20181113-GATK4/ANNOTATE_VCF/Ce330_annotated.vcf.gz |\
	bcftools query -i 'GT="alt"' -s N2 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk '$3 != "0/1" {print}' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_N2_GATK4_masked_no_het.txt
	rm CBbe.sqlite
done;

for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../input_files/Ce330_GATK4_STRELKA2_Intersect.vcf.gz |\
	bcftools query -i 'GT="alt"' -s N2 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_N2_GATK4.txt
	rm CBbe.sqlite
done;

for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../input_files/Ce330_GATK4_STRELKA2_Intersect.vcf.gz |\
	bcftools query -i 'GT="alt"' -s CB4851 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_CB4851_GATK4.txt
	rm CBbe.sqlite
done;


for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../20181113-GATK4/ANNOTATE_VCF/Ce330_annotated.vcf.gz |\
	bcftools query -i 'GT="alt"' -s CB4851 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_CB4851_GATK4_masked.txt
	rm CBbe.sqlite
done;

for chrom in I II III IV V X; do
	bcftools filter -i N_MISSING=0 -r $chrom ../20181113-GATK4/ANNOTATE_VCF/Ce330_annotated.vcf.gz |\
	bcftools query -i 'GT="alt"' -s CB4856 -f '[%CHROM\t%POS\t%GT\n]' |\
	awk '$3 != "0/1" {print}' |\
	awk -F '\t' 'BEGIN{printf("create table T(c text,p int); BEGIN TRANSACTION;\n");} {printf("insert into T(c,p) values(\"%s\",%d);\n",$1,$2);} END {printf("COMMIT; SELECT c,CAST(p/1E4  AS INT)*1E4, count(*) from T group by c,CAST(p/1E4  AS INT)*1E4 ; drop table T;");}' |\
	sqlite3 -separator ' ' CBbe.sqlite > $chrom.1E4_CB4856_GATK4_masked_no_het.txt
	rm CBbe.sqlite
done;
