library(tidyverse)
library(data.table)
library(ape)

orthologs <- fread('final_gene.emapper.seed_orthologs',skip = 5,fill =T)
colnames(orthologs)[1] <- 'qseqid'
orthologs <- orthologs %>% mutate(qseqid = str_extract(qseqid,'\\w+'))

gene_stat <- as.data.frame(table(orthologs$sseqid)) %>% 
  filter(Var1 != '') %>%
  arrange(desc(Freq))

dup_gene <- gene_stat[gene_stat$Freq > 1 ,]$Var1
uni_gene <- gene_stat[gene_stat$Freq == 1 ,]$Var1

orthologs_dup <- orthologs[orthologs$sseqid %in% dup_gene ,]
orthologs_dup <- orthologs_dup %>% group_by(sseqid) %>% arrange(sstart)
fwrite(orthologs_dup[,1:2],file = 'duplicated_gene.txt',sep = '\t')

orthologs_unique <- orthologs[orthologs$sseqid %in% unique_gene]
fwrite(orthologs_unique[,1:2],file = 'unique_gene.txt',sep = '\t')

## construct new fasta file-----
trinity_res <- read.FASTA('../trinity/Trinity.fasta')
names(trinity_res) <- str_extract(names(trinity_res),'\\w+')
unique_gene <- trinity_res[orthologs_unique$qseqid]
names(unique_gene) <- orthologs_unique$sseqid
write.FASTA(unique_gene,file = 'unique_gene.fa')
###
duplicate_gene <- trinity_res[orthologs_dup$qseqid]
names(duplicate_gene) <- orthologs_dup$sseqid
write.FASTA(duplicate_gene,file = 'duplicate_gene.fa')