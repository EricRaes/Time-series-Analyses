library(tidyverse)
library(phyloseq)
library(microbiome)
library(reshape2)
library(ggpubr)


#########################
# Import data
#######################
(Bedford<-read_csv2phyloseq(
  otu.file = "./Data_files/Bedford_Basin/ASV_table.csv",
  taxonomy.file = "./Data_files/Bedford_Basin/Taxonomy.csv",
  metadata.file = "./Data_files/Bedford_Basin/Metadata.csv",
  sep = ","))

## Remove taxa
Bedford_no_mito = subset_taxa(Bedford, !Kingdom=="NA")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Kingdom=="Eukaryota")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Phylum=="NA")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Order=="Chloroplast")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria")
Bedford_no_mito = subset_samples(Bedford_no_mito, Depth < 50) # remove deep samples





######################################
#Alpha diversity      
########################################
# use (or not) sample.size = 5000, to your read depth of interest
# 5000, 8000, 10000, 20000
Bedford_rar <-rarefy_even_depth(Bedford_no_mito, sample.size = 20000,
                                rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)

### for the rare microbiome, ASVs contributing less < 1% 
## if you rarefy to  5000 reads, then look at ASVs which contribute less than 50 reads. 50 is 1% of 5000 reads
(Bedford_rar= prune_taxa(taxa_sums(Bedford_rar) < 200, Bedford_rar))

# Compute metrics
#results <- estimate_richness(Bedford_rar, measures =c( 'Chao1', "Shannon")) # Phyloseq doesnt have Pielou
results <- alpha(Bedford_rar, index=c("chao1", "shannon", "pielou")) %>% # From microbiome package
  rownames_to_column("X") %>% 
  rename(Chao1=chao1, Shannon=diversity_shannon, Pielou=evenness_pielou)

metadata <- sample_data(Bedford_rar) %>% 
  data.frame() %>% 
  rownames_to_column("X")

METADATA <- merge(results,metadata, by="X", all=TRUE)

METADATA$Date_Merge <- as.Date(METADATA$Date_Merge)




###################################
#### Plotting alpha metrics over time
################################
my.cols_alphab <- c("green", "saddlebrown", "blue1", "slateblue",
                    "orange", "darkorchid", "red", 
                    "cyan", "Blue1", "black", 
                    "gold1", "violet", "goldenrod")
## Set factor levels
METADATA$Depth..m. = factor(METADATA$Depth..m., levels=c("1m", "5m", "10m"), 
                            labels=c("euphotic","euphotic", "euphotic"))
METADATA$Month = factor(METADATA$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))


METADATA$weekmonth <- paste(METADATA$Week, METADATA$Month, sep = "-")
data <- METADATA[c(2:4, 69)]
datlong = melt(data, id.vars = "weekmonth")
datlong[c('week', 'month')] <- str_split_fixed(datlong$weekmonth, '-', 2)
datlong$week <- as.numeric(datlong$week)


ggplot(datlong, aes(x= week, y= value, color= month)) +
  theme_bw()+geom_point(shape=19, size = 2) +
  facet_wrap(~ variable, scales = "free") +
  stat_smooth(method = "loess",  aes(color = "blue"), formula = y ~ x, se = TRUE) +
  scale_colour_manual(values=my.cols_alphab) +
  labs(y = "Alpha Diversity Measure") +
  theme(axis.title.x = element_text(size=20, vjust = 2),
        axis.title.y = element_text(size=20, vjust = 2),
        axis.text.y = element_text(size=20, vjust = 0.5, face="bold"),
        axis.text.x = element_text(size=20, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.title=element_text(size=13,face="bold"),
        legend.position="none",
        strip.text.x = element_text(size = 20),
        panel.spacing = unit(0.7, "lines"),
        panel.border = element_rect(colour = "black", linewidth=2),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.2, "cm"))

ggsave("output/Bedford_Basin/Chao_20000_rar.png", units = "cm" , height = 10, width = 20, dpi = 300)





###################################
#### Correlating alpha metrics chla, SST, and daylength
################################
## Set factor levels
METADATA$Depth..m. = factor(METADATA$Depth..m., levels=c("1m", "5m", "10m"), 
                            labels=c("euphotic","euphotic", "euphotic"))
METADATA$Month = factor(METADATA$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

# Chlorophyll.A  day_length Temperature
# Chao1 Shannon Pielou
## ggplot
ggplot(data=METADATA, aes(x=Temperature, y=Pielou, fill = "blue1")) +
  geom_point(aes(color="blue1"), shape = 21,size = 2,colour = "black") +
  scale_fill_manual(values="blue1") +
  scale_colour_manual(values="blue1") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=22, vjust = 0.5),
        axis.text.x = element_text(size=22, vjust = 0, angle = 0, hjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  stat_smooth(method=lm)+
  stat_cor(method = "pearson", aes(color = "blue1"), size=6, p.accuracy = 0.001, r.accuracy = 0.01)+
  guides(fill = "none", color = "none", linetype = "none", shape = "none")
##
#
#
#
#

ggsave("output/Bedford_Basin/pielou_Temperature.png")
