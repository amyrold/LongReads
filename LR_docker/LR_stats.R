#Final R Code for LongReads

# Folder Management ----
wk.dir <- getwd() 
# names of folders for output data (figures + data output)
folder.names <- c('1_raw_data', '2_filtered_data', '3_output', '4_blast')
# if folder with name "i" does not exist, create it.
for(i in 1:length(folder.names)){
  if(file.exists(folder.names[i]) == FALSE){
    dir.create(folder.names[i]) 
  }
}

# paths to the folders. The 'p.' indicates the variable is a path.
# make sure the variable names describe the folder.names
p_raw_data <- paste(wk.dir, "/", folder.names[1], "/", sep = "")
p_filt_data <- paste(wk.dir, "/", folder.names[2], "/", sep = "")
p_out <- paste(wk.dir, "/", folder.names[3], "/", sep = "")
p_blast<- paste(wk.dir, "/", folder.names[4], "/", sep = "")


#importing libraries
library(tidyverse)
library(ggplot2)
library(dplyr)


#Reading in the within and between files
within = read.csv(paste(p_out,"within.csv", sep = ""))
header = c("label", "ED")
colnames(within) <- header
within$label <- "Intravariation"

between = read.csv(paste(p_out,"between.csv", sep = ""))
header = c("label", "ED")
colnames(between) <- header
between$label <- "Intervariation"

#TEXT OUTPUT
  
sink(paste(p_out,"output.txt", sep = ""))

"Intervariation Distribution Summary"
#Intervariation distribution summary
summary(between$ED)  

cat("\n")
"Intravariation Distribution Summary"
#Intravariation distribution summary
summary(within$ED)

cat("\n")
"T-test of between (intergenomic) and within (intragenomic) edit distances"
#t-test 
t.test(between$ED, within$ED)

"Intervariation Breakdown of Shared and Unique Copies"
both <- rbind(within,between)
both %>% mutate(Type = case_when(ED > 0 ~ "Unique", ED == 0 ~ "Shared")) %>% group_by(label, Type) %>% summarize("count" = n()) %>% filter(label == "Intervariation") %>% mutate("total" = sum(count)) %>% mutate("percent" = (count/total)*100)
cat("\n")
"Intravariation Breakdown of Shared and Unique Copies"
both %>% mutate(Type = case_when(ED > 0 ~ "Unique", ED == 0 ~ "Shared")) %>% group_by(label, Type) %>% summarize("count" = n()) %>% filter(label == "Intravariation") %>% mutate("total" = sum(count)) %>% mutate("percent" = (count/total)*100)

sink()


#FIGURES

#Side by Side Boxplots
both <- rbind(within,between)

png(file=paste(p_out,"boxplots.png", sep = ""), res = 120, width=600, height=350)
ggplot(both, aes(x=label, y=ED, fill = label)) + geom_boxplot() + geom_violin(alpha = 0.5) + ggtitle('Boxplot of the Distribution of Edit Distances') 
dev.off()

#Distribution of copies
blast <- read.csv(paste(p_blast,"myresults.csv", sep = ""))
header <- c("accession", "pident", "qstart", "qend", "length", "strand", "evalue")
colnames(blast) <- header

png(file=paste(p_out,"copies_per_genome_dist.png", sep = ""), res = 120, width=600, height=350)
blast %>% filter(length > 1400) %>% group_by(accession) %>% summarize(num_copies = n()) %>%  ggplot(aes(x = num_copies)) + geom_histogram() + stat_bin(bins = 10) + ggtitle("Distribution of the number of 16S copies") + ylab("num_genomes")
dev.off()


#Shared and Unique
#unique is not 0, and shared is 0
#How many are shared and unique within intra and how many are shared and unique within inter

both <- rbind(within,between)
png(file=paste(p_out,"inter_intra_shared_unique.png", sep = ""), res = 120, width=600, height=350)
both %>% mutate(Type = case_when(ED > 0 ~ "Unique", ED == 0 ~ "Shared")) %>% group_by(label, Type) %>% summarize("count" = n()) %>% ggplot(aes(x = label, y = count, fill = Type)) + geom_bar(position = "fill", stat = "identity") + ylab("relative abundance") + ggtitle("Distribution of Shared and Unique Copies")
dev.off()
