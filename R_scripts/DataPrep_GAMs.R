library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R") 
#install_github("microbiome/microbiome") 
#devtools::install_github("vmikk/metagMisc")
library(qiime2R)
library(phyloseq)
library(breakaway)
library(microbiome)
library(metagMisc)
library(geosphere)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
#     Import data
################################################################################

##########
# A. Bedford Basin
##########
# Read file locations
otu_mat <-"../Data_files/Bedford_Basin/ASV_table.csv"
tax_mat <- "../Data_files/Bedford_Basin/Taxonomy.csv"
samples_df <- "../Data_files/Bedford_Basin/Metadata.csv"

# You can check for if file exists
file.exists(c(tax_mat, samples_df , otu_mat)) 

# Load 3 csv files into 1 phyloseq object (all in 1 Bedford Basion dataset)
(Bedford<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ","))

# Filter in the same way as in R_Scripts_Time_Series.R
(Bedford = subset_taxa(Bedford, !Kingdom=="NA"))
(Bedford = subset_taxa(Bedford, !Kingdom=="Eukaryota"))
(Bedford = subset_taxa(Bedford, !Phylum=="NA"))
(Bedford = subset_taxa(Bedford, !Order=="Chloroplast"))
(Bedford = subset_taxa(Bedford, !Family=="Mitochondria"))
(Bedford = subset_samples(Bedford, Depth < 50)) # remove deep samples

# Additionally:
(Bedford = prune_samples(sample_sums(Bedford) >= 5000, Bedford)) 

sample_data(Bedford)$Month = factor(sample_data(Bedford)$Month, levels=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sep","Oct", "Nov", "Dec"))
sample_data(Bedford)$Year. = factor(sample_data(Bedford)$Year., levels=c("Year_2014", "Year_2015", "Year_2016", "Year_2017"))

Bedford_break <- breakaway(Bedford, cutoff = 10)  # Computes a model to predict true richness

## Set factor levels

Bedford_break_summary <- Bedford_break %>% summary
names(Bedford_break_summary)[5] <- "SampleID"
META<- read.csv("../Data_files/Bedford_Basin/Metadata.csv", header=TRUE)

#######merge with Metadata
Bedford_break_summary<-merge(META, Bedford_break_summary, by="SampleID")
Bedford_break_summary$Depth..m. = factor(Bedford_break_summary$Depth..m., levels=c("1m", "5m", "10m"),
                                         labels=c("euphotic","euphotic", "euphotic"))

write_csv(Bedford_break_summary, "../Data_files/Bedford_Richness_Meta.csv")


##########
# B. L4, English Channel
##########
# Read file locations
otu_mat <-"../Data_files/L4_Engl_Channel/ASV_table.csv"
tax_mat <- "../Data_files/L4_Engl_Channel/Taxonomy.csv"
samples_df <- "../Data_files/L4_Engl_Channel/Metadata.csv"

# Load 3 csv files into 1 phyloseq object (all in 1 Bedford Basion dataset)
(L4_EC <-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ","))

L4_EC_break <- breakaway(L4_EC)

L4_EC_break <- L4_EC_break %>% summary
names(L4_EC_break)[5] <- "SampleID"

# import metadata
L4_EC_meta <- read.csv("../Data_files/L4_Engl_Channel/Metadata.csv", header=TRUE)

# Merge metadata with breakaway stats and select only important columns
L4_EC_meta_break <- merge(L4_EC_meta, L4_EC_break, by="SampleID") %>% 
  as_tibble()

# Cast "Character" column as "date"
L4_EC_meta_break$datetime <- as.Date(L4_EC_meta_break$Date)

# Sort dataset by time
L4_EC_meta_break <- arrange(L4_EC_meta_break, Date)

# Check
view(L4_EC_meta_break)

write_csv(L4_EC_meta_break, "../Data_files/L4_Richness_Meta.csv")

##########
# C. Blaynes Bay Microbial Observatory (BBMO)
##########

# loading Metadata
metadata <-  read_csv("../Data_files/BBMO/Metadata.csv") %>%
  column_to_rownames(var = "SampleID")     # make sample ID the row names
metadata$Year <- as.integer(metadata$Year)  # Cast data as integer for later ploting
metadata$weeknum <- as.integer(metadata$weeknum)
# loading ASV table
ASV_table <- read_csv("../Data_files/BBMO/ASV_table.csv") %>%
  column_to_rownames('ASV_NAME')      # Make ASV IDs row names

# Create 1 phyloseq object with 2 components
(BBMO <-  phyloseq(otu_table(ASV_table, taxa_are_rows = TRUE), sample_data(metadata)))

BBMO_break <- breakaway(BBMO, cutoff = 15)

# Make a table relating breakaway estimates to SampleIDs
BBMO_break <- BBMO_break %>% summary
names(BBMO_break)[5] <- "SampleID"

# import metadata (again?)
BBMO_meta <- read.csv("../Data_files/BBMO/Metadata.csv", header=TRUE) %>% 
  

# Merge metadata with breakaway stats and select only important columns
BBMO_meta_break <- merge(BBMO_meta, BBMO_break, by="SampleID") %>% 
  as_tibble()  

# Creating date column by arbitrarily assigning the 1st of each month as the date of sampling (unavailable info)
BBMO_meta_break$datetime <-  as.Date(with(BBMO_meta_break, paste(Year, monthnum, "01", sep="-")), "%Y-%m-%d")

# Sort dataset by time
BBMO_meta_break <- arrange(BBMO_meta_break, datetime)

# Check
view(BBMO_meta_break)

write_csv(BBMO_meta_break, "../Data_files/BBMO_Richness_Meta.csv")

##########
# D. San Pedro Ocean Time Series (SPOTS)
##########

# Loading Metadata
metadata <-  read_csv("../Data_files/SPOTS/Metadata_5m.csv") %>%
  column_to_rownames(var = "SampleID") 
  
rownames(metadata) <- str_replace(rownames(metadata), "SPOTS", "SPOT")
# Loading ASV table
ASV_table <- read_csv("../Data_files/SPOTS/ASV_table_5m.csv") %>%
  column_to_rownames('ASV_id')      # Make ASV IDs row names

#tax <- read_csv("../Data_files/SPOTS/taxonomy.csv") %>%
#  column_to_rownames('ASV_id')      # Make ASV IDs row names

# Create 1 phyloseq object with 2 components
(SPOTS <- phyloseq(otu_table(ASV_table, taxa_are_rows = TRUE), 
                   #tax_table(tax),
                   sample_data(metadata))) 

SPOTS_break <- breakaway(SPOTS)

# Make a table relating breakaway estimates to SampleIDs
SPOTS_break <- SPOTS_break %>% summary
names(SPOTS_break)[5] <- "SampleID"

# import metadata (again?)
SPOTS_meta <- read.csv("../Data_files/SPOTS/Metadata_5m.csv", header=TRUE)

SPOTS_meta$SampleID <- str_replace(SPOTS_meta$SampleID, "SPOTS", "SPOT")

# Merge metadata with breakaway stats and select only important columns
SPOTS_meta_break <- merge(SPOTS_meta, SPOTS_break, by="SampleID") %>% 
  as_tibble() 

# Creating date column
SPOTS_meta_break$datetime <-  as.Date(with(SPOTS_meta_break, paste(year, month_num, day, sep="-")), "%Y-%m-%d")

# Sort dataset by time
SPOTS_meta_break <- arrange(SPOTS_meta_break, datetime)

# manually calculate daylength because this isn't in metadata supplied from github
# Lat from manuscript 33º30’
SPOTS_meta_break$daylength <- geosphere::daylength(33.3, SPOTS_meta_break$datetime)

ggplot(SPOTS_meta_break, aes(daylength, estimate))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

# Check
view(SPOTS_meta_break)

write_csv(SPOTS_meta_break, "../Data_files/SPOTS_Richness_Meta.csv")


##########
# E. Fram Strait
##########

# Loading Metadata
metadata <-  read_csv("../Data_files/Fram_Strait/Metadata.csv") %>%
  column_to_rownames(var = "SampleID")    # make sample ID the row names
# Loading ASV table
ASV_table <- read_csv("../Data_files/Fram_Strait/ASV_table.csv")[1:5906,] %>% 
  column_to_rownames('ASV')      # Make ASV IDs row names

# Create 1 phyloseq object from 3 components
(FRAM <- phyloseq(otu_table(ASV_table, taxa_are_rows = TRUE), sample_data(metadata))) 

FRAM_break <- breakaway(FRAM, cutoff = 15)

# Make a table relating breakaway estimates to SampleIDs
FRAM_break <- FRAM_break %>% summary
names(FRAM_break)[5] <- "SampleID"

# import metadata (again?)
FRAM_meta <- read.csv("../Data_files/Fram_Strait/Metadata.csv", header=TRUE)

# Merge metadata with breakaway stats
FRAM_meta_break <- merge(FRAM_meta, FRAM_break, by="SampleID") %>% 
  as_tibble() 

# Cast "Character" column as "date"
FRAM_meta_break$datetime <- as.Date(FRAM_meta_break$date)

# Sort dataset by time
FRAM_meta_break <- arrange(FRAM_meta_break, date)

# Check
view(FRAM_meta_break)

write_csv(FRAM_meta_break, "../Data_files/FRAM_Richness_Meta.csv")


##########
# F.G.H. Australian sites
##########
# Loading Metadata
metadata <-  read_csv("../Data_files/Australia/contextual_META.csv") %>%
  column_to_rownames(var = "SampleID")    # make sample ID the row names
# Loading ASV table
ASV_table <- read_csv("../Data_files/Australia/MAI_ROT_YON_ASV.csv") %>% 
  column_to_rownames('OTU')      # Make ASV IDs row names

# Create 1 phyloseq object from 3 components
(Australia <- phyloseq(otu_table(ASV_table, taxa_are_rows = TRUE), sample_data(metadata))) 

# Maria Island
(maria <- subset_samples(Australia, Station=='Maria Island'))
(maria <- filter_taxa(maria, function(x) sum(x) > 0 , TRUE))

# Yongala
(yongala <- subset_samples(Australia, Station=='Yongala'))
(yongala <- filter_taxa(yongala, function(x) sum(x) > 0 , TRUE))

# Rottnest Island
(rottnest <- subset_samples(Australia, Station=='Rottnest Island'))
(rottnest <- filter_taxa(rottnest, function(x) sum(x) > 0 , TRUE))

maria_break <- breakaway(maria, cutoff = 15)
yongala_break <- breakaway(yongala, cutoff = 15)
rottnest_break <- breakaway(rottnest, cutoff = 15)

##########
# F. Maria Island
##########
# Make a table relating breakaway estimates to SampleIDs
maria_break <- maria_break %>% summary
names(maria_break)[5] <- "SampleID"

# import metadata as data frame
maria_meta <-  read_csv("../Data_files/Australia/contextual_META.csv") %>% 
  filter(Station=="Maria Island")

# Merge metadata with breakaway stats and select only important columns
maria_meta_break <- merge(maria_meta, maria_break, by="SampleID") %>% 
  as_tibble() 

# Creating date column
maria_meta_break$datetime <-  as.Date(with(maria_meta_break, paste(Year, Month...20, Day, sep="-")), "%Y-%m-%d")

# Sort dataset by time
maria_meta_break <- arrange(maria_meta_break, datetime)

# Check
view(maria_meta_break)

write_csv(maria_meta_break, "../Data_files/Maria_Richness_Meta.csv")

##########
# G. Yongala
##########
# Make a table relating breakaway estimates to SampleIDs
yongala_break <- yongala_break %>% summary
names(yongala_break)[5] <- "SampleID"

# import metadata as data frame
yongala_meta <-  read_csv("../Data_files/Australia/contextual_META.csv") %>% 
  filter(Station=="Yongala")

# Merge metadata with breakaway stats and select only important columns
yongala_meta_break <- merge(yongala_meta, yongala_break, by="SampleID") %>% 
  as_tibble() 

# Creating date column
yongala_meta_break$datetime <-  as.Date(with(yongala_meta_break, paste(Year, Month...20, Day, sep="-")), "%Y-%m-%d")

# Sort dataset by time
yongala_meta_break <- arrange(yongala_meta_break, datetime)

# Check
view(yongala_meta_break)

write_csv(yongala_meta_break, "../Data_files/Yongala_Richness_Meta.csv")

##########
# H. Rottnest Island
##########
# Make a table relating breakaway estimates to SampleIDs
rottnest_break <- rottnest_break %>% summary
names(rottnest_break)[5] <- "SampleID"

# import metadata as data frame
rottnest_meta <-  read_csv("../Data_files/Australia/contextual_META.csv") %>% 
  filter(Station=="Rottnest Island")

# Merge metadata with breakaway stats and select only important columns
rottnest_meta_break <- merge(rottnest_meta, rottnest_break, by="SampleID") %>% 
  as_tibble() 

# Creating date column
rottnest_meta_break$datetime <-  as.Date(with(rottnest_meta_break, paste(Year, Month...20, Day, sep="-")), "%Y-%m-%d")

# Sort dataset by time
rottnest_meta_break <- arrange(rottnest_meta_break, datetime)

# Check
view(rottnest_meta_break)

write_csv(rottnest_meta_break, "../Data_files/Rottnest_Richness_Meta.csv")

