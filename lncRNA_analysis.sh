#! /bin/bash
# @Ajeet Singh
# singh.ajeet@nipgr.ac.in

fastq_dump=/home/user1/softwares/sratoolkit.2.10.2-centos_linux64/bin/fastq-dump
fastp=/home/user1/softwares/fastp
hisat2=/home/user1/softwares/hisat2-2.1.0/hisat2
stringtie=/home/user1/softwares/stringtie-2.1.1.Linux_x86_64/stringtie

hisat_index=/home/user1/ajeet/gossipium_barb/genome/gb_index
gtf=/home/user1/ajeet/gossipium_barb/genome/GCA_008761655.1_Gossypium_barbadense_v1.1_genomic.gtf
SRRacc=$1

$fastq_dump ${SRRacc} -e 30 --split-3

name=`basename "${SRRacc}"|sed 's/.sra//g'|tr '\n' ' '|sed 's/ $//g'`

num_file=`ls -1|wc -l|tr '\n' ' '|sed 's/ $//g'`

if [[ $num_file=>3 ]]
then
${fastp} -i ${SRRacc}\_1.fastq -I ${SRRacc}\_2.fastq -o ${SRRacc}\_filtered_1.fastq -O ${SRRacc}\_filtered_2.fastq -w 16 --detect_adapter_for_pe
${hisat2} --novel-splicesite-outfile novel_splicesite.txt --summary-file hisat_summary.txt --dta -x ${hisat_index} -1 $1\_filtered_1.fastq -2 $1\_filtered_2.fastq -p 35 -S $name\.sam

elif [[ $num_file==2 ]]
then
${fastp} -i $1\.fastq -o $1\_filtered.fastq -w 16
${hisat2} --novel-splicesite-outfile novel_splicesite.txt --summary-file hisat_summary.txt --dta -x ${hisat_index} -U $1\_filtered.fastq -p 35 -S $name\.sam
fi

map=`tail -n 1 hisat_summary.txt |awk -F'%| ' '{print $1}'|tr '\n' ' '|sed 's/ $//g'`
if [[ "$map" > "50" ]]
then
samtools sort -@ 40 $name\.sam -o $name\.bam
${stringtie} $name\.bam -G ${gtf} -o $name\.gtf -p 40 -l $name -A $name\-gene.tab

elif [[ "$map" < "50" || "$map" == "50" ]]
then
printf "\n============\nMapping Quality of data is bad!!!!\n=============\n"
fi

mkdir -p reads && mv *sra* reads
rm -f reads/*fastq

#find ./ -name '*.gtf' >ALL_GTF.txt
#/home/user1/softwares/stringtie-2.1.1.Linux_x86_64/stringtie --merge -o merged.gtf -m 200 -p 60 ALL_GTF.txt
#sed -n '3,$'p merged.gtf|grep -v exon|awk -F'\t|"' '{print $1,$4,$5,$12,"1",$7}' OFS='\t' >transcript.bed
#Install cpc2, python2, virtualenv -p /usr/bin/python2 -p env, source env/bin/activate
#source ~/softwares/CPC2-beta/bin/env/bin/activate
#(env) [user1@localhost gossipium_barb]$ python ~/softwares/CPC2-beta/bin/CPC2.py -i merge.fasta -o cpc2-out

#source ~/softwares/CPC2-beta/bin/env/bin/activate ..same env activate
#(env) [user1@localhost gossipium_barb]$ python ~/softwares/CNCI/CNCI.py -f merge.fasta -o cnci-out.txt -p 40 -m pl
#/home/user1/softwares/CNCI
#CNCI classification were completely done!
#384813.609675 second for 220291 transcript's computation.

# Running plncpro:
# python prediction.py -p plncpro-result-file -i merge.fasta -m models/dicot.model -o gossipiyum-out -d lib/blastdb/db -t 15 -r
# [user1@localhost zFINAL-set]$ grep -vf pfam-coding <(comm -12 <(sort CNCI-lnc) <(sort cpc2-lnc)|comm -12 - <(sort plncpro-lnc))|wc -l
# 86914
# [user1@localhost zFINAL-set]$ comm -12 <(sort CNCI-lnc) <(sort cpc2-lnc)|comm -12 - <(sort plncpro-lnc)|wc -l
# 86933

#Pfam search
#sed '/>/s/:/_/g' ../gh_merge.fasta >change_header_gh_merge.fasta
#/home/user1/softwares/EMBOSS-6.6.0/bin/transeq change_header_gh_merge.fasta gh_merge.pep -frame 6
#hmmsearch --cpu 50 -o gh_out --tblout gh_tblout --domtblout gh_domtblout --pfamtblout gh_pfamtblout /home/user1/ajeet/gossipium_barb/Pfam/Pfam-A.hmm ../gh_merge.pep

#comm -12 <(cut -f1 cnci-noncode|sort) <(cut -f1 cpc2-noncode|sort)|comm -12 <(cut -f1 plncpro-noncode|sort) - >aa
#comm -31 <(sort Pfam-hits-id) <(sort aa) >lncrna-id
#stringtie -e -p 18 -G stringtie_merged.gtf -o re-estimation-stringtie/$i/$i.gtf $i/$i.bam
