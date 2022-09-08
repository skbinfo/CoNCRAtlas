# CoNCRAtlas :An atlas of cotton non-coding RNAs
https://nipgr.ac.in/CoNCRAtlas
<a href="https://doi.org/10.5281/zenodo.7057078"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7057078.svg" alt="DOI"></a>
<img src="http://14.139.61.8/CoNCRAtlas/images/lncRNA_pipeline.png" width="700" height="700">

__lncRNA_pipe.sh__ is a pipeline for identification, expression and classification of lncRNAs from multiple RNA-Seq datasets.
Input required genome, annotation file (GTF), and multiple RNA-Seq datasets.
It gives lncRNAs - fasta, GTF, gemonic classification, gene count, and gene expression files in output.

```
git clone https://github.com/skbinfo/CoNCRAtlas
cd CoNCRAtlas
chmod a+x lncRNA_pipe.sh
```
Usage:
```
lncRNA_pipe.sh -h

   -x  <hisat2-build indexes>
   -i  <fastq file; single-end or _1.fastq if paired >
	     for multiple single-end files or forward pair of paired-end use it multiple time
		 -i first.fq -i second.fq .. or -i first_1.fq -i second_1.fq ..
   -I  <fastq file; _2.fastq paired >
	     for multiple reverse pair of paired-end files use it multiple time
		 -I first_2.fq -I second_2.fq ..
   -g  <genome file>
   -a  <annotation file; GTF>
   -o  <outdir>

  Miscellanious:
   -t  <int> number of threads [default:1]
   -r  <yes/no> remove intermediate files [yes or no; default: no]

```
Usage example:
For Single-end reads
```
lncRNA_pipe.sh -i <foo.fq> -i <faa.fq> -g <genome> -a <genome GTF> -o <out-dir>
```
For Paired-end reads
```
lncRNA_pipe.sh -i <foo_1.fq> -i <faa_1.fq> -I <foo_2.fq> -I <faa_2.fq> -g <genome> -a <genome GTF> -o <out-dir>
```

