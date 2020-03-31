library(DESeq2)
library(tidyverse)
library(magrittr)
options(digits=2)

#readcount <- read.csv("/Users/chentao/Desktop/Build/ReadCount.csv",
#                      header = TRUE, stringsAsFactors = FALSE)
#group <- list(control=c("Dev_3M_1", "Dev_3M_2", "Dev_3M_3"),
#              treat=c("Dev_1M_1", "Dev_1M_2", "Dev_1M_3"))
DESeq_s <- function(object, group) {
  #build data matrix
  cts = object %>% 
    dplyr::select(2:ncol(object)) %>%
    `rownames<-`(object[,1]) %>% 
    as.matrix
  #build group data.frame
#  coldata = group %>% as.data.frame %>% 
#    gather("control", "treat", key="condition", value="sample",
#           factor_key=TRUE) %>% 
#    `rownames<-`(.[[2]]) %>% 
#    select("condition")
  sum_group <- summary(group)
  coldata = data.frame(condition=rep(rownames(sum_group), sum_group[,"Length"]),
                       sample=unlist(group))
  rownames(coldata) <- coldata[,"sample"]
  
  #build DEseq object
  rdds <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ condition)
  #start analysis
  dds <- DESeq(rdds)
  res <- results(dds)
  
  #get DESeq result
  normalized1 <- res@listData %>% as.data.frame %>%
    mutate(GeneID=res@rownames) %>%
    dplyr::select(GeneID, baseMean, log2FoldChange, pvalue, padj) %>% na.omit
  #get normalzed value
  res1 <- counts(dds, normalized=TRUE)
  normalized2 <- res1 %>% as.data.frame(., row.names = FALSE) %>%
    mutate(GeneID=rownames(res1))
  
  #build result
  rrsult <- list(readcount=object,
                 diff=normalized1,
                 normalized=normalized2)
  return(rrsult)
}