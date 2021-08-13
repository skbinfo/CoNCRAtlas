#gb_SRA_sheet (SRA data downloaded from NCBI and filtered manually for wild type samples only)
awk -F"\t" 'BEGIN {OFS="\t"}  $1 == "" { $i = "NA" }  {print}' gb_SRA_sheet > gh_pheno.csv (To fill blank cells with NA)

#Get IDs of lncRNAs predicted by FEELnc
awk '{print $3}' gb_lncrna_classification_isBest.txt > FEELnc_IDs

#expression (FPKM step remaining)
sh chk_para_lncRNA.sh FEELnc_IDs gh_pheno.csv
awk -F\\t '{if($2!=0)print}' gb_lncrna_SRR_expression.csv >gb_lncrna_SRR_expression_filter.csv
awk -F"\t" '!seen[$1, $5]++' gb_lncrna_SRR_expression_filter.csv >gb_lncrna_SRR_expression_filter_keep_uni.csv
awk '{print $1}' gb_express.csv |sort -u > gb_ids
awk 'FNR==NR{a[$1];next}($1 in a){print}' gb_ids gb_lncrna_SRR_expression_filter_keep_uni.csv > gb_express.csv
     # add FPKM cutoff step
sed -f dbIDs_gb.sed gb_express.csv > gb_express_new.csv (To add DBIDs)

#final gtf
grep -wf gb_ids merged.gtf

#final fasta
xargs samtools faidx merge_edit.fasta < gb_ids >> gb_lncRNA.fa

#Classification
grep -wf gb_ids gb_lncrna_classification_isBest.txt > gb_lncrna_Best.csv
cut --complement -f1 gb_lncrna_Best.csv >gb_lncrna_classfied.csv
cat <(head -1 gb_lncrna_classfied.csv|sed 's/^/DBID\t/g') <(sed -n '2,$'p gb_lncrna_classfied.csv|sort|awk -F\\t -v a="0001" '{printf "%s%05d%s\n","CoLNCGB",a,"\t"$0;a+=1}' OFS="\t" ) >gb_lnc_classification.csv


#Analysis of tools results

 #CNCI
grep "noncoding" CNCI.index > CNCI_nc.index
awk '{print $1}' CNCI_nc.index > CNCI_nc_IDs.index
cat gb_ids | while read i;do  t=$(echo "$i"|awk -F: '{print $1}'); grep -w "$t" CNCI_nc.index ; done >> CNCI_final.index

 #cpc2
grep "noncoding" cpc2-out.txt > cpc2_nc.txt
awk '{print $1}' cpc2_nc.txt > cpc2_nc_IDs.txt
cat gb_ids | while read i;do  t=$(echo "$i"|awk -F: '{print $1}'); grep -w "$t" cpc2_nc.txt ; done >> cpc2_final.txt

 #plncpro
awk '{if($2 == 0){print}}' plncpro-result-file > plncpro-nc-gb
awk '{print $1}' plncpro-nc-gb > plncpro-nc-gb-IDs
cat fpkm_more_0.1| while read i;do  t=$(echo "$i"|awk -F: '{print $1}'); grep -w "$t" plncpro-nc-gb ; done >> plncpro_final.txt

 #pfam
grep "MSTRG" gh_domtblout |awk -F'[\t;_]' '{ print $1,$3,$4}'|awk '{ print $1 "::" $2 ":" $3 }' > pfam_gb_final (have duplicates)
sort pfam_gb_final | awk '!seen[$0]++' > pfam_IDs.txt (no duplicates)


 #Finding Common IDs of CNCI, cpc2, plncpro & pfam
sort CNCI_nc_IDs.index > sorted_CNCI_IDs.txt
sort cpc2_nc_IDs.txt > sorted_cpc2_IDs.txt
sort plncpro-nc-gb-IDs > sorted_plncpro_IDS.txt
comm -12 sorted_cpc2_IDs.txt sorted_CNCI_IDs.txt 
comm -12 sorted_cpc2_IDs.txt sorted_CNCI_IDs.txt > common_in_2
comm -12 common_in_2 sorted_plncpro_IDS.txt 
comm -12 common_in_2 sorted_plncpro_IDS.txt > common_in_3
comm -23 common_in_3 pfam_IDs.txt > total_com_IDs.txt

#fasta file Common IDs of CNCI, cpc2, plncpro & pfam
xargs samtools faidx gh_merge.fasta < total_com_IDs_gh.txt >> fasta_gh.fas

#removing strands with unknown type
grep -v "(\.)" total_com_IDs_gh.txt < strands_IDs_gh

#GTF of Common IDs of CNCI, cpc2, plncpro & pfam
step1:
cat merged.gtf|while read i; do id=$(echo "$i"|grep -o "transcript_id.*"|sed 's/;.*//g;s/transcript_id "//g;s/"//g') ; add=$(echo "$i"|awk '{print "::"$1":"$4"-"$5"("$7")"}'); printf "$id$add\t$i\n" ; done >> gb_merged.gtf

step2:
cat strands_IDs.txt| while read i;do  t=$(echo "$i"|awk -F: '{print $1}'); grep -w "$t" gb_merged.gtf; done >>gb_final.gtf
