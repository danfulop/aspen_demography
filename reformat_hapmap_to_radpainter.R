# reformat hapmap to radpainter

library(dplyr)
library(tidyr)

hapmap <- read.table("HapMap.hmp.8base_filtered.mod.renamed.txt", 
                     header=T, check.names=F, stringsAsFactors=F)
hapmap <- hapmap[-c(3:11)] # drop columns I don't need
hapmap %>% names
hapmap[1:6,1:6] %>% summary
hapmap[1:20,1:4] 

switch_base <- function(x, alleles) {
  ifelse(x=="N", "", 
         ifelse(x=="A" | x=="T" | x=="G" | x=="C", x, alleles))
}

switch_base(test$alleles, test$`1KS_Keizer_F_OR`) # works vectorized!

rad <- hapmap

rad[3:ncol(rad)] <- apply(rad[3:ncol(rad)], 2, switch_base, rad$alleles)

rad[1:20, 1:4]

write.table(rad[3:ncol(rad)], file="rad_input.tsv", quote=F, sep="\t", row.names=F)
