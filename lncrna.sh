#Get transcripts IDs
grep ">" gh_merge.fasta |sed 's/>//g'

#CNCI
grep "noncoding" CNCI.index > CNCI_nc.index
awk '{print $1}' CNCI_nc.index > CNCI_nc_IDs.index

#cpc2
grep "noncoding" cpc2-out.txt > cpc2_nc.txt
awk '{print $1}' cpc2_nc.txt > cpc2_nc_IDs.txt

#plncpro
awk '{if($2 == 0){print}}' plncpro-result-file > plncpro-nc-gb
awk '{print $1}' plncpro-nc-gb > plncpro-nc-gb-IDs

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

#fasta file of total common IDs
xargs samtools faidx gh_merge.fasta < total_com_IDs_gh.txt >> final_fasta_gh.fas

#removing strands with unknown type
grep -v "(\.)" total_com_IDs_gh.txt < strands_IDs_gh

#GTF from fasta file
step1:
cat merged.gtf|while read i; do id=$(echo "$i"|grep -o "transcript_id.*"|sed 's/;.*//g;s/transcript_id "//g;s/"//g') ; add=$(echo "$i"|awk '{print "::"$1":"$4"-"$5"("$7")"}'); printf "$id$add\t$i\n" ; done >> gb_merged.gtf

step2:
cat strands_IDs.txt| while read i;do  t=$(echo "$i"|awk -F: '{print $1}'); grep -w "$t" gb_merged.gtf; done >>gb_final.gtf
