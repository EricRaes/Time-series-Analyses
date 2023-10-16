################
# Linear regression
###############
library(tidyverse)
library(phyloseq)
library(microbiome)
library(relaimpo)

### FRAM
# Metadata
metadata <-  read_csv("Data_files/Fram_Strait/Metadata.csv") %>%
  column_to_rownames(var = "SampleID") %>% # make sample ID the row names
  mutate(Week = strftime(date, format = "%V")) 
#metadata$Year <- as.integer(metadata$Year)  # Cast data as integer for later ploting
metadata$Week <- as.integer(metadata$Week)
# ASV table (Sample X Taxon, read counts)
ASV_table <- read_csv("Data_files/Fram_Strait/ASV_table.csv") %>%
  column_to_rownames('ASV')      # Make ASV IDs row names
# Making phyloseq object
(Fram <-  merge_phyloseq(otu_table(ASV_table, taxa_are_rows = TRUE), 
                         sample_data(metadata))) 



######################################
#Alpha diversity and Correlations     #change sample.size = 5000, to your read depth of interest
########################################

Fram_rar <-rarefy_even_depth(Fram, sample.size = 10842,
                                rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)

#results <- estimate_richness(Bedford_rar, measures =c( 'Chao1', "Shannon")) # Phyloseq doesnt have Pielou
results <- alpha(Fram_rar, index=c("chao1", "shannon", "pielou")) %>% # From microbiome package
  rownames_to_column("X") %>% 
  rename(Chao1=chao1, Shannon=diversity_shannon, Pielou=evenness_pielou)




d <- sample_data(Fram_rar) %>% 
  as.data.frame() %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var="X")


d.all<-merge(results,d, by="X", all=TRUE)
d.all <- d.all %>% rename(Day = daylight)
d.all <- d.all %>% rename(Temp = temp) 
d.all <- d.all %>% rename(Chl = chl_sens) 


fit_chao1 <- lm(Chao1 ~  Temp + Chl + Day, data = d.all)
fit_Shannon <- lm(Shannon ~ Temp + Chl + Day , data = d.all)

summary(fit_chao1)
step(fit_chao1) 
residuals(fit_chao1)

###
# diagnostic plots
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(fit_Shannon)

# Calculate Relative Importance for Each Predictor
calc.relimp(fit_Shannon,
            type=c("lmg","last","first","pratt"),
            rela=TRUE)

# Bootstrap Measures of Relative Importance (1000 samples)
boot <- boot.relimp(fit_chao1, 
                    b = 1000, 
                    type = c("lmg", "last", "first", "pratt"), 
                    rank = TRUE,
                    diff = TRUE, 
                    rela = TRUE)

booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result

