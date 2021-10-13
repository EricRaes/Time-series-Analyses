if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
devtools::install_github("adw96/breakaway", force = T)
devtools::install_github("adw96/DivNet")
install.packages("remotes")
remotes::install_github("YuLab-SMU/ggtreeExtra")
BiocManager::install("ComplexHeatmap")
library(tidyverse)
library(qiime2R)
library(ggplot2)
library(microbiome); packageVersion("microbiome")
library(breakaway)
library(DivNet)
library(magrittr)
library(phyloseq)
library(ggpubr)

setwd("~/Documents/DivNet/Eric/DivNet")

##import data

otu_mat <-"ASV_BB_Table.csv"
tax_mat <- "Taxonomy_species_Silva138.1_BB.csv"
samples_df <- "METADATA_niskin_CTD_R_3.csv"

file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

carbom<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")
##########clean up
carbom_no_mito = subset_taxa(carbom, !Order=="Chloroplast" )
carbom_no_mito = subset_taxa(carbom_no_mito, !Family=="Mitochondria" )
carbom_no_mito = subset_samples(carbom_no_mito, !Month== "XXX")
carbom_no_mito = subset_samples(carbom_no_mito, !Depth..m.== "60m")
##############
carbom_no_mito= prune_samples(sample_sums(carbom_no_mito) >= 5000, carbom_no_mito) #clean up your data and remove samples with low sequencing depth; in this case only keep samples which have a sequencing depth >= 5000 reads
carbom_no_mito= subset_samples(carbom_no_mito, sample_names(carbom_no_mito) != "BB14-22C") # Remove sample with poor model fit in breakaway
#####

#This is more computationally intensive (typical 1-min runtime)
samp_rich <- sample_richness(carbom_no_mito)

freq <- carbom_no_mito %>% otu_table %>% 
  build_frequency_count_tables

# breakaway function for cutoff value
cutoff_wrap <- function(my_data, requested = NA) {
  
  iss <- my_data$index
  fis <- my_data$frequency
  length_fis <- length(fis)
  
  breaks <- which(iss[-1] - iss[-length_fis] > 1)
  cutoff <- ifelse(is.na(breaks[1]), length_fis, breaks[1])
  
  if (!is.na(requested)) {
    # if the requested cutoff is lte cutoff, use cutoff, o/w, cutoff
    if (requested <= cutoff) {
      cutoff <- requested # ok
    } else {
      warning("ignoring requested cutoff; it's too low")
      cutoff <- cutoff
    }
  }
  
  return(cutoff)
}

est_rich <- breakaway(carbom_no_mito, cutoff = 15)

#Some samples have extreme uncertainties based on the frequency counts
plot(est_rich, carbom_no_mito, color = "Year")

summary(est_rich) %>% as_tibble() 

#Model selection was inspected for specific samples (e.g., BB14-22C with large error bars and sequence depth = 12,387) 
plot(est_rich$'BB14-22C')
#Here, a transformed weighted linear regression model(tWLRLM), instead of the Kemp models described in Willis et al (2015) was fit in a sample with high singleton count that inflates predicted frequency of unobserved taxa 
#Experimenting with function 'breakaway_nof1' to remove singletons fixes high error but misfits numerous other samples. However, put in perspective poor model fit to BB14-22C leaves overall patterns unaffected.

##### summary and dataframe it

break_summary <- est_rich %>% summary
names(break_summary)[5] <- "SampleID"
META<- read.csv("METADATA_niskin_CTD_R_3.csv", header=TRUE)

#######merge with Metadata
brokenaway<-merge(META, break_summary, by="SampleID")
brokenaway$Depth..m. = factor(brokenaway$Depth..m., levels=c("1m", "5m", "10m"),
                            labels=c("euphotic","euphotic", "euphotic"))

#Create new variable for seasons to conduct hypothesis testing
brokenaway$Season <- brokenaway$Month
#Winter
brokenaway$Season[brokenaway$Season =="Dec"]<- "Winter"
brokenaway$Season[brokenaway$Season =="Jan"]<- "Winter"
brokenaway$Season[brokenaway$Season =="Feb"]<- "Winter"
brokenaway$Season[brokenaway$Season =="Mar"]<- "Winter"
#Spring
brokenaway$Season[brokenaway$Season =="Apr"]<- "Spring"
brokenaway$Season[brokenaway$Season =="May"]<- "Spring"
#Summer
brokenaway$Season[brokenaway$Season =="Jun"]<- "Summer"
brokenaway$Season[brokenaway$Season =="Jul"]<- "Summer"
brokenaway$Season[brokenaway$Season =="Aug"]<- "Summer"
#Autumn
brokenaway$Season[brokenaway$Season =="Sep"]<- "Autumn"
brokenaway$Season[brokenaway$Season =="Oct"]<- "Autumn"
brokenaway$Season[brokenaway$Season =="Nov"]<- "Autumn"

##plot it
p<-ggplot(data=brokenaway,aes(x=Week, y=estimate, colour=Year.))
p <- p + theme_bw(24)+
  geom_point(shape = 1,size = 2)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, size=1)+
#  stat_regline_equation(aes(label = ..rr.label..), size=8)+
  geom_jitter(width = 0, height = 0)+
  scale_color_manual(values=c("blue1", "orange","#06f5fb",
                              "black" ,"red", "#a7cc7b",
                              "#cba9e5", "violet", "#808a8a",
                              "orange", "#ff249d" ))+
  #facet_grid(~Depth)+
  theme(axis.title.x = element_text(size=20, vjust = 0.3),
        axis.title.y = element_text(size=20, vjust = 0.9),
        axis.text.y = element_text(size=20, vjust = 0.9),
        axis.text.x = element_text(size=20, vjust = 0.3, angle = 0, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))+
  theme(panel.grid.major = element_blank())+
  stat_smooth(method = "loess", aes(color = Depth..m.),
              formula = y ~ x, se = TRUE)+
  ggtitle("Breakaway")+
  ylab("Richness")+
  ylim(0,900)

ggsave("Breakaway.pdf", p, height = 5 , width = 5* asp, dpi = 300)

#Hypothesis testing for richness
my_x <- model.matrix(~Season, data = brokenaway)[, -1]
bm <- breakaway::betta(chats = betta.df.min$Estimate, ses = betta.df.min$Error, X = my_x)

bt <- betta(chats = summary(est_rich)$estimate,
            ses = summary(est_rich)$error,
            X= my_x)
bt$table

########### DivNet ########### 
#lp <- carbom_no_mito %>% tax_glom("Genus")

div<-divnet(carbom_no_mito, 
                      base = "893ddb0ea7f30889b34ae3cfd79137eb",
                      variance = "parametric",
                      ncores = 4)

#DivNet has high comp. demands (use divnet rust version with large datasets)

##### summary and dataframe it
divnet_summary <- divnet_phylum$shannon %>% summary
names(divnet_summary)[5] <- "SampleID"

#######merge with Metadata
DivNeted<-merge(META, divnet_summary, by="SampleID")
DivNeted$Depth..m. = factor(DivNeted$Depth..m., levels=c("1m", "5m", "10m"),
                            labels=c("euphotic","euphotic", "euphotic"))
##plot it
div <- ggplot(data = DivNeted, aes(x= Week, y=estimate, color= Year.))+
  geom_point(shape = 1,size = 2.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, size=1)+
  geom_jitter(width = 0, height = 0)+
  scale_color_manual(values=c("blue1", "orange","#06f5fb",
                              "black" ,"red", "#a7cc7b",
                              "#cba9e5", "violet", "#808a8a",
                              "orange", "#ff249d" ))+
  #facet_grid(~Depth)+
  theme_bw(24)+
  theme(axis.title.x = element_text(size=20, vjust = 0.3),
        axis.title.y = element_text(size=20, vjust = 0.9),
        axis.text.y = element_text(size=20, vjust = 0.9),
        axis.text.x = element_text(size=20, vjust = 0.3, angle = 0, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))+
  theme(panel.grid.major = element_blank())+
  stat_smooth(method = "loess", aes(color = Depth..m.),
              formula = y ~ x, se = TRUE)+
  ggtitle("DivNet")+
  ylab("Shannon Estimate")

ggsave("Shannon_Proteobacteria_whole.pdf", div, height = 5 , width = 5* asp, dpi = 300)

#Hypothesis testinf with divnet estimates
my_x <- model.matrix(~Season, data = DivNeted)
bt <- betta(chats = divnet_summary$estimate,
            ses = divnet_summary$error,
            X= my_x)
bt$table
###########
