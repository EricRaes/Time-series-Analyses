library(phyloseq)
library(microbiome)
library(tidyverse)
library(reshape2)
library(vegan)


######## Import BB data and make phyloseq
(Bedford<-read_csv2phyloseq(
  otu.file = "./Data_files/Bedford_Basin/ASV_table.csv",
  taxonomy.file = "./Data_files/Bedford_Basin/Taxonomy.csv",
  metadata.file = "./Data_files/Bedford_Basin/Metadata.csv",
  sep = ","))


######## Filter data
(Bedford_no_mito = subset_taxa(Bedford, !Kingdom=="NA"))
(Bedford_no_mito = subset_taxa(Bedford_no_mito, !Kingdom=="Eukaryota"))
(Bedford_no_mito = subset_taxa(Bedford_no_mito, !Phylum=="NA"))
(Bedford_no_mito = subset_taxa(Bedford_no_mito, !Order=="Chloroplast"))
(Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria"))
# Select depths
(Bedford_1 = subset_samples(Bedford_no_mito, Depth == 1)) # remove deep samples 
(Bedford_1_5_10 = subset_samples(Bedford_no_mito, Depth < 50)) # remove deep samples 
# Sum depths at each date
(Bedford_1_5_10 = merge_samples(Bedford_1_5_10, "Date_Merge", fun = sum))


######## Compute count of days for vegdist to compute euclidean distance in time
# Extract sample IDs and associated dates - order by date
meta_1 <- sample_data(Bedford_1) %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  select(c("Sample", "Date_Merge")) %>%
  mutate(Date_Merge = as.Date(Date_Merge)) %>% 
  arrange(Date_Merge)
# Generate a variable (Days) which is a count of days since first sample (day 0)
startdate <- meta_1$Date_Merge[1] # start date (sample 1)
meta_1 <- mutate(meta_1, Days = as.numeric(difftime(meta_1$Date_Merge, startdate, units = "days"))) %>% 
  select(c("Sample", "Days")) # Grab Sample ID and day count


######## Grab ASV table and transpose so taxa are columns
asv_table_1 <- otu_table(Bedford_1) %>% 
  as.data.frame() %>% 
  t()


######## Order rows of both ASV table and metadata 
######## This is done alphabetically on Sample IDs, not necessarily chronologically
# Order ASV table
asv_table_1 <- asv_table_1[order(row.names(asv_table_1)), ] 
# Order metadata the same way
meta_1 <- meta_1[order(meta_1$Sample), ] 



##########################################################
# Aitchison distances Time lag analysis
##########################################################
# Plot time Lag regression using Aitchison and Bray-Curtis
aitchison <- data.frame(aitchison = as.numeric(vegdist(asv_table_1, method = "aitchison", pseudocount=1)), time.lag = as.numeric(vegdist(meta_1$Days, method = "euclidean")))
bray.curtis <- data.frame(bray.curtis = as.numeric(vegdist(asv_table_1, method = "bray")), time.lag = as.numeric(vegdist(meta_1$Days, method = "euclidean")))

aitchison$time.lag <- round(aitchison$time.lag, digits = -1)
bray.curtis$time.lag <- round(bray.curtis$time.lag, digits = -1)

bray.curtis %>% 
  #mutate(Lag = sqrt(interval)) %>%
  ggplot(aes(x = time.lag, y = bray.curtis)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  #geom_point() +
  #stat_smooth(method = "lm", se = F, linewidth = 1) + # fit a regression
  #stat_smooth(method = "loess", se = F, size = 1, col = "firebrick3") + # fit a curve
  xlab("Time lag (days)") + ylab("Aitchison distance") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_line(colour = "grey75", linewidth = 0.5),
        plot.margin=unit(c(0.2,1,0.2,0.2),"cm"))

ggsave("output/Bedford_Basin/bb_5m_aitchison_tl.png", units = "cm" , height = 15, width = 25, dpi = 300)







#########################
# Ordination PCA CLR and Euclidian distance
#######################
######### Set colours
my.cols <- c("darkorchid","orange", "blue1",
             "green" ,"black",
             "cyan", "red", "saddlebrown",
             "goldenrod", "violet",
             "gold1", "slateblue", "chartreuse")

########## Set factor levels
sample_data(Bedford_no_mito)$Month = factor(sample_data(Bedford_no_mito)$Month, levels=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sep","Oct", "Nov", "Dec"))
sample_data(Bedford_no_mito)$Year. = factor(sample_data(Bedford_no_mito)$Year., levels=c("Year_2014", "Year_2015", "Year_2016", "Year_2017"))

##### Compute aitchison distances to RDA
set.seed(66)
Bedford_clr <- microbiome::transform(Bedford_no_mito, "clr")   
out.pcoa.logt <- ordinate(Bedford_clr, method = "RDA", distance = "euclidean")
evals <- out.pcoa.logt$CA$eig

### Plot
p3<-plot_ordination(Bedford_clr, out.pcoa.logt, type = "Sample", color = "Month") 
#p3$layers <- p3$layers[-1]
p3 + #ggtitle("PCA Bedford 16S rRNA ASV")+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.7))+
  theme(legend.key=element_blank())+
  geom_point(size=3)+
  scale_colour_manual(values=my.cols)
