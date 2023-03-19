# Termite_transcriptome
This is the bioinfomatics analysis pipeline for termite_transcriptome paper

## *De novo* assembly RNA-seq reads
>zcat ./clean_data/*_1.clean.fq.gz |bgzip -@ 40 > ./trinity/all_left.fq.gz \
zcat ./clean_data/*_2.clean.fq.gz |bgzip -@ 40 > ./trinity/all_right.fq.gz


>Trinity --seqType fq --left ./trinity/all_left.fq.gz 
--right ./trinity/all_right.fq.gz 
--min_contig_length 300 --output trinity 
--CPU 40 --max_memory 200G --include_supertranscripts 

*De novo assembly* with [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)

>TransDecoder.LongOrfs -t Trinity.fasta 

Predict ORF with [TransDecoder](https://github.com/TransDecoder/TransDecoder)

>cd-hit -i longest_orfs.fa -o unique_pep.fa -c 0.9 -n 5 -T 40 -d 100

Reduce sequence redundancy by [cdhit](https://github.com/weizhongli/cdhit) clustering threshold is the 90% identity.

>./emapper.py -m diamond --cpu 40 --report_orthologs -i unique_pep.fa -o final_gene

Annotation gene by [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)

Run Get_unique_gene.R to get unique_gene.fa and duplicate_gene.fa.
>python merge.py duplicate_gene.fa > duplicate_gene.fa \
>cat unique_gene.fa duplicate_gene.fa > all_gene.fa



>./salmon-1.9.0_linux_x86_64/bin/salmon index -t all_gene.fa -i termites_index -k 21

Construct salmon index
>snakemake -s salmon_quant.smk -c 40

Quantification gene expression by [salmon](https://github.com/COMBINE-lab/salmon).

You can contact the authors of the manuscript (kongyunhui1@gmail.com) with any questions or requests for code, which will be addressed as soon as possible.
