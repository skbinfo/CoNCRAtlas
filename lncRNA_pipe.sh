#! /bin/bash
# @Ajeet Singh
# singh.ajeet@nipgr.ac.in

# This script identify the confident lncRNAs from RNA-Seq datasets

#command args 
while getopts 'x:i:I:g:a:o:t:r:h' c
do
	case $c in
		x) index="$OPTARG";;
		i) reads+=("$OPTARG");;
		I) READS+=("$OPTARG");;
		g) genome="$OPTARG";;
		a) annotation="$OPTARG";;
		o) outdir="$OPTARG";;
		t) threads="$OPTARG";;
		r) remove_i="$OPTARG";;
		h) echo ""
		   echo "	This pipeline gives the confident lncRNAs from multiple RNA-Seq datasets. It gives lncRNAs fasta file, annotation file,"
		   echo "	genomic classification, and lncRNA expression for each sample"
		   echo "	"
		   echo "   -x  <hisat2-build indexes>"
		   echo "   -i  <fastq file; single-end or _1.fastq if paired >"
		   echo "	     for multiple single-end files or forward pair of paired-end use it multiple time"
		   echo "		 -i first.fq -i second.fq .. or -i first_1.fq -i second_1.fq .."
		   echo "   -I  <fastq file; _2.fastq paired >"
		   echo "	     for multiple reverse pair of paired-end files use it multiple time"
		   echo "		 -I first_2.fq -I second_2.fq .."
		   echo	"   -g  <genome file>"
		   echo	"   -a  <annotation file; GTF>"
		   echo "   -o  <outdir>"
		   echo ""
		   echo "  Miscellanious:"
		   echo "   -t  <int> number of threads [default:1]"
		   echo "   -r  <yes/no> remove intermediate files yes or no"
		   echo ""
		   exit 1;;
	   [?]) printf "\n      usage: ./lncRNA_pipe.sh -x [hisat_indexes] -i [fastq reads] -o [output dir] \n\n"
		   exit 1;;	
	esac
done
shift $((OPTIND-1))

#check command & tools

fastp=/home/user1/softwares/fastp
hisat2=/home/user1/softwares/hisat2-2.1.0/hisat2
eexon=/home/user1/softwares/hisat2-2.1.0/extract_exons.py
ess=/home/user1/softwares/hisat2-2.1.0/extract_splice_sites.py
hisat2_build=/home/user1/softwares/hisat2-2.1.0/hisat2-build
stringtie=/home/user1/softwares/stringtie-2.1.1.Linux_x86_64/stringtie
bedtools=bedtools
transeq=/home/user1/softwares/EMBOSS-6.6.0/bin/transeq
plncpro_dir=/home/user1/softwares/plncpro
classify=/home/user1/softwares/FEELnc/FEELnc_classifier.pl
featureCounts=/home/user1/softwares/subread-2.0.1-Linux-x86_64/bin/featureCounts

if ! command -v ${hisat2} &> /dev/null
then
        echo "${hisat2} not found"
        exit 1;
elif ! command -v ${fastp} &> /dev/null
then
        echo "${fastp} not found or not in path"
        exit 1;
elif ! command -v ${stringtie} &> /dev/null
then
	echo "${stringtie} not found or not in path"
	exit 1;
elif ! command -v ${bedtools} &> /dev/null
then
	echo "${bedtools} not found or not in path"
	exit 1;
elif ! command -v ${transeq} &> /dev/null
then
	echo "${transeq} not found or not in path"
	exit 1;
elif ! command -v ${classify} &> /dev/null
then
	echo "${classify} not found or not in path"
	exit 1;
elif ! command -v ${featureCounts} &> /dev/null
then
        echo "${featureCounts} not found or not in path"
        exit 1;
fi

# virtualenv -p /usr/bin/python2 -p myenv
# Install cpc2, CNCI
activate_myenv () {
source /home/user1/softwares/myenv/bin/activate
}

lib_exp () {
export PERL5LIB=$PERL5LIB:/home/user1/ajeet/gossipium_barb/zFINAL-set/FEELnc/lib/
}

hisat_index_dir=/home/user1/ajeet/gossipium_barb/genome
hisat_index=/home/user1/ajeet/gossipium_barb/genome/gb_index
gtf=/home/user1/ajeet/gossipium_barb/genome/GCA_008761655.1_Gossypium_barbadense_v1.1_genomic.gtf

if [ `ls ${hisat_index}* 2>/dev/null | wc -l ` -eq 0 ]
then
	python2 ${eexon} ${gtf} >${hisat_index_dir}/exon.exon
	python2 ${ess} ${gtf} >${hisat_index_dir}/splice.ss
	${hisat2_build} --ss ${hisat_index_dir}/splice.ss --exon ${hisat_index_dir}/exon.exon ${genome} ${hisat_index_dir}/gb_index
fi

if [ ${reads} ] && [ ${READS} ]
then
        for inda in "${!reads[@]}"; do
                for indb in "${!READS[@]}"; do
                        if [ ${inda} -eq ${indb} ]
                        then
                                echo "${reads[$inda]} ${READS[$indb]}"
				bn=$(basename ${reads[$inda]} _1.fastq)
				${fastp} -i ${reads[$inda]} -I ${READS[$indb]} -o ${outdir}/${bn}_filtered_1.fastq -O ${outdir}/${bn}_filtered_2.fastq -w ${threads} --detect_adapter_for_pe
				rem="${outdir}/${bn}_filtered_1.fastq ${outdir}/${bn}_filtered_2.fastq ${outdir}/fastp.json ${outdir}/fastp.html"
				${hisat2} --summary-file ${outdir}/${bn}_hisat_summary.txt --dta -x ${hisat_index} -1 ${outdir}/${bn}_filtered_1.fastq -2 ${outdir}/${bn}_filtered_2.fastq -p ${threads} -S ${outdir}/${bn}.sam
				rem+=" ${outdir}/${bn}_hisat_summary.txt ${outdir}/${bn}.sam"
				map=`tail -n 1 ${outdir}/${bn}_hisat_summary.txt |awk -F'%| ' '{print $1}'|tr '\n' ' '|sed 's/ $//g'`
				if [[ "$map" > "50" ]]
				then
					samtools sort -@ ${threads} ${outdir}/${bn}.sam -o ${outdir}/${bn}.bam
					rem+=" ${outdir}/${bn}.bam"
					${stringtie} ${outdir}/${bn}.bam -G ${gtf} -o ${outdir}/${bn}.gtf -p ${threads} -l $bn -A ${outdir}/${bn}-gene.tab
					rem+=" ${outdir}/${bn}.gtf ${outdir}/${bn}-gene.tab"
				elif [[ "$map" < "50" || "$map" == "50" ]]
				then
					printf "\n============\nMapping Quality of data is bad!!!!\n=============\n"
				fi
			fi
                done
        done
elif [ $reads ] && [ -z $READS ]
then
        for sseq in "${reads[@]}"; do
                echo "$sseq"
				bn=$(basename $sseq .fastq)
				${fastp} -i $sseq -o ${outdir}/${bn}_filtered.fastq -w ${threads}
				rem="${outdir}/${bn}_filtered.fastq"
				${hisat2} --summary-file ${outdir}/${bn}_hisat_summary.txt --dta -x ${hisat_index} -U ${outdir}/${bn}_filtered.fastq \
					-p ${threads} -S ${outdir}/${bn}.sam
				rem+=" ${outdir}/${bn}_hisat_summary.txt ${outdir}/${bn}.sam"
				map=`tail -n 1 ${outdir}/${bn}_hisat_summary.txt |awk -F'%| ' '{print $1}'|tr '\n' ' '|sed 's/ $//g'`
				if [[ "$map" > "50" ]]
				then
					samtools sort -@ ${threads} ${outdir}/${bn}.sam -o ${outdir}/${bn}.bam
					rem+=" ${outdir}/${bn}.bam"
					${stringtie} ${outdir}/${bn}.bam -G ${gtf} -o ${outdir}/${bn}.gtf -p ${threads} -l $bn -A ${outdir}/${bn}-gene.tab
					rem+=" ${outdir}/${bn}.gtf ${outdir}/${bn}-gene.tab"
				elif [[ "$map" < "50" || "$map" == "50" ]]
				then
					printf "\n============\nMapping Quality of data is bad!!!!\n=============\n"
				fi
        done
fi

find ${outdir} -name '*.gtf' > ${outdir}/ALL_GTF.txt
rem+=" ${outdir}/ALL_GTF.txt"
${stringtie} --merge -o ${outdir}/merged.gtf -m 200 -p ${threads} ${outdir}/ALL_GTF.txt
rem+=" ${outdir}/merged.gtf"
sed -n '3,$'p merged.gtf|grep -v exon|awk -F'\t|"' '{if($7!=".")print $1,$4,$5,$12,"1",$7}' OFS='\t' >${outdir}/merged.bed
rem+=" ${outdir}/merged.bed"
${bedtools} getfasta -fi ${genome} -bed ${outdir}/merged.bed -s -name -fo ${outdir}/merged.fasta
rem+=" ${outdir}/merged.fasta"

activate_myenv
python ~/softwares/CPC2-beta/bin/CPC2.py -i ${outdir}/merged.fasta -o ${outdir}/cpc2_out.txt
python ~/softwares/CNCI/CNCI.py -f ${outdir}/merged.fasta -o ${outdir}/cnci_out.txt -p $threads -m pl
deactivate
rem+=" ${outdir}/cpc2_out.txt ${outdir}/cnci_out.txt"

grep "noncoding" ${outdir}/cpc2_out.txt | awk '{print $1}' | sort > ${outdir}/sorted_cpc2.IDs
grep "noncoding" ${outdir}/cnci_out.txt | awk '{print $1}' | sort > ${outdir}/sorted_cnci.IDs
rem+=" ${outdir}/sorted_cpc2.IDs ${outdir}/sorted_cnci.IDs"

# Running plncpro:
#download the swissprot database and make blast database of that by makblastdb to lib/blastdb/db in plncpro
python ${plncpro_dir}/prediction.py -p ${outdir}/plncpro_out.txt -i ${outdir}/merged.fasta \
	-m ${plncpro_dir}/models/dicot.model -o ${outdir} -d ${plncpro_dir}/lib/blastdb/db -t $threads -r
rem+=" ${outdir}/plncpro_out.txt"

awk '{if($2 == 0)print}' ${outdir}/plncpro_out.txt | awk '{print $1}' |sort > ${outdir}/sorted_plncpro.IDs
rem+=" ${outdir}/sorted_plncpro.IDs"

#Pfam search
sed '/>/s/:/_/g' ${outdir}/merged.fasta >${outdir}/merged_change_header.fasta
${transeq} ${outdir}/merged_change_header.fasta ${outdir}/merged_change_header.pep -frame 6
rem+=" ${outdir}/merged_change_header.fasta ${outdir}/merged_change_header.pep"

#Download Pfam and make search the translated peptides for domains
#wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
${hmmsearch} --cpu $threads -o ${outdir} --tblout pfam_tblout --domtblout pfam_domtblout --pfamtblout pfam_pfamtblout \
	/home/user1/Pfam/Pfam-A.hmm ${outdir}/merged_change_header.pep
rem+=" ${outdir}/pfam_tblout ${outdir}/pfam_domtblout ${outdir}/pfam_pfamtblout"

grep "MSTRG" ${outdir}/pfam_domtblout |awk -F'[\t;_]' '{ print $1,$3,$4}'|awk '{ print $1 "::" $2 ":" $3 }' \
	| sort | awk '!seen[$0]++' > ${outdir}/sorted_pfam.IDs
rem+=" ${outdir}/sorted_pfam.IDs"

comm -12 ${outdir}/sorted_cpc2.IDs ${outdir}/sorted_cnci.IDs > ${outdir}/common_cpc2_cnci
comm -12 ${outdir}/common_cpc2_cnci ${outdir}/sorted_plncpro.IDs > ${outdir}/common_cpc2_cnci_plncpro
comm -23 ${outdir}/common_cpc2_cnci_plncpro ${outdir}/sorted_pfam.IDs > ${outdir}/only_noncode.IDs
rem+=" ${outdir}/common_cpc2_cnci ${outdir}/common_cpc2_cnci_plncpro ${outdir}/only_noncode.IDs"

sed 's/\./\\./g' ${outdir}/only_noncode.IDs > ${outdir}/only_noncode_for_gtf.IDs
rem+=" ${outdir}/only_noncode_for_gtf.IDs"
cat ${outdir}/only_noncode_for_gtf.IDs|while read i
do
	grep -w "$i" ${outdir}/merged.gtf
done | awk '!seen[$0]++' > ${outdir}/lncRNA.gtf

lib_exp
${classify} -i ${outdir}/lncRNA.gtf -a ${gtf} > ${outdir}/lncRNA_classes.txt

if [ ${reads} ] && [ ${READS} ]
then
        for inda in "${!reads[@]}"; do
                for indb in "${!READS[@]}"; do
                        if [ ${inda} -eq ${indb} ]
                        then
                                echo "${reads[$inda]} ${READS[$indb]}"
				bn=$(basename ${reads[$inda]} _1.fastq)
				${stringtie} -e -p $threads -G ${outdir}/lncRNA.gtf -o ${outdir}/${bn}.re.gtf ${outdir}/${bn}.bam
                        fi
                done
        done
elif [ $reads ] && [ -z $READS ]
then
        for sseq in "${reads[@]}"; do
                echo "$sseq"
		bn=$(basename $sseq .fastq)
		${stringtie} -e -p $threads -G ${outdir}/lncRNA.gtf -o ${outdir}/${bn}.re.gtf ${outdir}/${bn}.bam
        done
fi

all_bam=$(find ${outdir} -name '*.bam'|tr '\n' ' ')
$featureCounts -a ${gtf} -o ${outdir}/lncRNA_gene_count_exon.out ${all_bam} -T $threads -p
$featureCounts -a ${gtf} -o ${outdir}/lncRNA_gene_count_exonMO.out ${all_bam} -T $threads -p -M -O --fraction

if [ "$remove_i" == "[Nn][Oo]" ]
then
	:
elif [ "$remove_i" == "[Yy][Ee][Ss]" ]
then
	rm -rf $rem
fi
