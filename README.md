# CoNCRAtlas :An atlas of cotton non-coding RNAs
__lncRNA_pipe.sh__ is a pipeline for identification, expression and classification of lncRNAs from multiple RNA-Seq datasets.
Input required genome, annotation file (GTF), and multiple RNA-Seq datasets.
It gives lncRNAs - fasta, GTF, gemonic classification, gene count, and gene expression files in output.
```
git clone https://github.com/skbinfo/CoNCRAtlas
cd CoNCRAtlas
chmod a+x lncRNA_pipe.sh
```
Usage example:
```
lncRNA_pipe.sh -i <foo.fq> -i <faa.fq> -g <genome> -a <genome GTF>
```


<img src="http://14.139.61.8/CoNCRAtlas/images/lncRNA_pipeline.png" width="700" height="700">
