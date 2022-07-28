library(jsonlite)
library(purrr)
library(RCurl)
library(tidyverse)
library(AnnotationForge)
update_kegg <- function(json = "ko00001.json") {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  
  kegg <- fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <-
          kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <-
          str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <-
          rbind(pathway2name,
                tibble(Pathway = pathway_id, Name = pathway_name))
        
        kos_info <-
          kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        
        kos <- str_match(kos_info, "K[0-9]*")[, 1]
        
        ko2pathway <-
          rbind(ko2pathway, tibble(
            Ko = kos,
            Pathway = rep(pathway_id, length(kos))
          ))
      }
    }
  }
  save(pathway2name, ko2pathway, file = "kegg_info.RData")
}
update_kegg(json = "ko00001.json")
##update kegg information
##https://www.genome.jp/kegg-bin/get_htext?ko00001 download ko00001.json first
##code from http://www.genek.tv/course/225/task/4861/show

orthologs_ano <- data.table::fread('final_gene.emapper.annotations',skip = 4,fill =T) %>%
  distinct(seed_ortholog,GOs,KEGG_ko) %>%
  filter(seed_ortholog != '') %>% 
  as.data.frame()

gene_info <- orthologs_ano %>% 
  dplyr::select(GID = seed_ortholog, GENENAME = seed_ortholog) %>% 
  distinct(GID,GENENAME)  

gene2go_info <- orthologs_ano %>% 
  dplyr::select(GID = seed_ortholog, GOs) %>% 
  distinct(GID,GOs) %>%
  filter(GID != '',GOs != '-')
gene2go  <- map_dfr(1:nrow(gene2go_info),function(x,gene_info){
  go <-  str_split(gene_info[x,2],',')[[1]]
  df <- data.frame(GID = rep(gene_info[x,][[1]],length(go)),
                   GO  =go,
                   EVIDENCE = rep("IEA", length(go)))
  return(df)
},gene_info = gene2go_info)
##Construct Gene to Go


load('kegg_info.RData')
gene2Ko_info <- orthologs_ano %>% 
  dplyr::select(GID = seed_ortholog, Ko = KEGG_ko) %>%
  distinct(GID,Ko) %>%
  filter(GID != '',Ko != '-')
gene2ko  <- map_dfr(1:nrow(gene2Ko_info),function(x,gene_info){
  ko <-  str_split(gene_info[x,2],',')[[1]]
  df <- data.frame(GID = rep(gene_info[x,][[1]],length(ko)),
                   Ko  =ko)
  return(df)
},gene_info = gene2Ko_info)
gene2ko$Ko <- str_replace(gene2ko$Ko,"ko:","")

gene2pathway <- gene2ko %>% 
  left_join(ko2pathway, by = "Ko") %>%
  dplyr::select(GID, Pathway) %>% 
  distinct(GID, Pathway) %>%
  filter(!is.na(Pathway))

# tax_id = "60588"
# genus = "Odontotermes" 
# species = "formosanus"

makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway= gene2pathway,
               version="0.0.1",
               maintainer = "kongyunhui <kongyunhui1@gmail.com>",
               author = "kongyunhui <kongyunhui1@gmail.com>", 
               outputDir = ".",  
               tax_id= '60588',  
               genus= 'Odontotermes',  
               species= 'formosanus', 
               goTable="go")
##Run in shell
##tar zcvf ./org.Oformosanus.eg.db.tar.gz ./org.Oformosanus.eg.db
install.packages("org.Oformosanus.eg.db.tar.gz", repos = NULL, type = "source")