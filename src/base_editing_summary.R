library("data.table")
library("ggplot2")
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)

sample_name = args[1]
base_editing_file = args[2]


tbl = fread(base_editing_file)

editingSummary = tbl %>% group_by(guide) %>% count(editType)
editingSummary = editingSummary %>% group_by(guide) %>% mutate(countT= sum(n)) %>%
  group_by(guide, add=TRUE) %>%
  mutate(per=paste0(round(100*n/countT,2),'%'))

write.table(editingSummary,file=paste(sample_name,"editing_summary.txt",sep="_"),sep="\t",quote=F,row.names = F)

editingPositions = tbl %>% group_by(guide) %>% count(editPositions)
editingPositions = editingPositions %>% group_by(guide) %>% mutate(countT= sum(n)) %>%
  group_by(guide, add=TRUE) %>%
  mutate(per=paste0(round(100*n/countT,2),'%'))
  
  
write.table(editingPositions,file=paste(sample_name,"editing_positions.txt",sep="_"),sep="\t",quote=F,row.names = F)