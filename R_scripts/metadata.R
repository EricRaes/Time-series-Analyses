library(plyr)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(GGally)
library(ggsci)
library(geosphere)
library(microbiome)

### BBMO
# Metadata
metadata <-  read_csv("Data_files/BBMO/Metadata.csv") %>%
  column_to_rownames(var = "SampleID") %>%    # make sample ID the row names
  mutate(ENV_CHL_total=log(ENV_CHL_total))
metadata$Year <- as.integer(metadata$Year)  # Cast data as integer for later ploting
metadata$weeknum <- as.integer(metadata$weeknum)
# ASV table (Sample X Taxon, read counts)
ASV_table <- read_csv("Data_files/BBMO/ASV_table.csv") %>%
  column_to_rownames('ASV_NAME')      # Make ASV IDs row names
# Making phyloseq object
(BBMO <-  merge_phyloseq(otu_table(ASV_table, taxa_are_rows = TRUE), 
                         sample_data(metadata))) 

# Rarefying 
BBMO_rar <-rarefy_even_depth(BBMO, sample.size = 17500,
                             rngseed = 123, replace = TRUE, 
                             trimOTUs = TRUE, verbose=TRUE) 

#results <- estimate_richness(BBMO_rar, measures =c( 'Chao1', "Shannon")) # Phyloseq doesnt have Pielou
results <- alpha(BBMO_rar, index=c("chao1", "shannon", "pielou")) %>% # From microbiome package
  rownames_to_column("X") %>% 
  dplyr::rename(Chao1=chao1, Shannon=diversity_shannon, Pielou=evenness_pielou)

# Extract metadata from rarefied dataset
metadata <- sample_data(BBMO_rar) %>% 
  data.frame() %>% 
  rownames_to_column("X")
# Merge div indices with metadata
BBMO_div_metadata <- merge(results,metadata, by="X", all=TRUE)

# Here is what I really care for
ofinterest <- BBMO_div_metadata[c('X',
                                  "Chao1",
                                  "Shannon",
                                  "Pielou",
                                  'ENV_Temp',
                                  'ENV_CHL_total',
                                  'ENV_PO4',
                                  'ENV_NH4',
                                  'ENV_NO2',
                                  'ENV_NO3',
                                  'ENV_SI',
                                  'ENV_Day_length_Hours_light')]

# correlogram and more
ggpairs(ofinterest[2:11], cardinality_threshold = NULL)






#################################################
# all bellow is code to make metadata plots with sites combined
#################################################

################
# Loading data
################
BB <- read_csv("Data_files/Bedford_Basin/Metadata.csv") %>% 
  dplyr::rename(Station=Voyage, dl=day_length, chla=Chlorophyll.A, sst=Temperature) %>% 
  dplyr::filter(Depth..m. %in% c("1m", "5m", "10m")) %>% 
  select(c(44, 6, 3, 19, 24, 20, 41)) #%>% 
  #mutate(chla=log(chla))
BB <- mutate(BB, NOx=BB$Nitrate + BB$Nitrite)%>% 
  select(c(1:5, 8))

BBMO <- read_csv("Data_files/BBMO/Metadata.csv")[c(21,16, 5, 2, 8, 9)] %>% 
  dplyr::rename(Week=weeknum, dl=ENV_Day_length_Hours_light, chla=ENV_CHL_total, sst=ENV_Temp) #%>% 
  #mutate(chla=log(chla))
BBMO$Station <- "Blanes Bay"
BBMO <- mutate(BBMO, NOx=BBMO$ENV_NO2 + BBMO$ENV_NO3) %>% 
  select(c(1:4, 7:8))

Fram <- read_csv("Data_files/Fram_Strait/Metadata.csv") %>% 
  mutate(Week = strftime(date, format = "%V")) %>% 
  mutate(Station="Fram Strait") %>% 
  dplyr::rename(dl=daylight, chla=chl_sens, sst=temp, NOx=NO3_NO2) %>% 
  select(c(Week, Station, dl, chla, sst, NOx)) #%>% 
  #mutate(chla=log(chla))

L4 <- read_csv("Data_files/L4_Engl_Channel/Metadata.csv")[c(6, 13, 17, 3, 14)] %>% 
  dplyr::rename(dl=day, chla="Chlorophyll A (ug/L)", sst="Temperature (C)", NOx="NO2 + NO3 (umol L-1)") #%>% 
  #mutate(chla=log(chla))
L4$Station <- "English Channel"

SPOT <- read_csv("Data_files/SPOTS/Metadata_5m.csv") %>% 
  dplyr::rename(Week=week_num) %>% 
  mutate(Station="San Pedro", 
         latitude=33.3,
         NOx=NA) #%>% 
  #mutate(chla=log(chla))
SPOT$date <- as.Date(with(SPOT, paste(year, month_num, day, sep="-")), "%Y-%m-%d") 
SPOT$dl <- daylength(SPOT$latitude, SPOT$date)
SPOT <- SPOT[c(8:10,12,14, 16)]

################
# northern hemisphere
################
northern <- rbind(SPOT, L4, Fram, BBMO, BB) %>% 
  drop_na(Station) %>% 
  mutate(Week=as.integer(Week))

N_colors <- c("#56B4E9", "#006600", "#990099", "#000066", "#E69F00")

################
# southern hemisphere
################
southern <- read_csv("Data_files/Australia/contextual_META.csv")[c(8, 12, 14, 74, 88, 36, 37)]%>% 
  dplyr::rename(dl=day_length, chla=CPHL_A, sst=CTDTemperature) %>% 
  drop_na(Station) #%>% 
  #mutate(chla=log(chla))
southern <- rowwise(southern) %>% 
  mutate(NOx = sum(c(Nitrate_umol_L, Nitrite_umol_L), na.rm = TRUE)) %>% 
  select(c(1:5, 8))


S_colors <- c("#33FF33", "#9B26B6" ,"#E10600")

################
# Temp
################
ggplot(northern, aes(x=Week, y= sst, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 2) +
  stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = "Temperature (°C)") + 
  scale_color_manual(values = N_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
  )#legend.position="none")
ggsave("output/Nouthern_leg.png", units = "cm" , height = 10, width = 10, dpi = 300)


ggplot(southern, aes(x=Week, y= sst, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 2) +
  stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = "Temperature (°C)") + 
  scale_color_manual(values = S_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position="none")

ggsave("output/Southern_temp.png", units = "cm" , height = 10, width = 10, dpi = 300)

################
# chla
################
dev.new()
ggplot(northern, aes(x=Week, y= chla, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 2) +
  stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = substitute(paste("Chlorophyll - ", italic("a"), " (μg L"^1, ")"))) + 
  scale_color_manual(values = N_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position="none")
ggsave("output/Nouthern_chla.png", units = "cm" , height = 10, width = 10, dpi = 300)


ggplot(southern, aes(x=Week, y= chla, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 2) +
  stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = substitute(paste("Chlorophyll - ", italic("a"), " (μg L"^1, ")"))) + 
  scale_color_manual(values = S_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position="none")

ggsave("output/Southern_chla.png", units = "cm" , height = 10, width = 10, dpi = 300)

################
# day length
################
ggplot(northern, aes(x=Week, y= dl, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 2) +
  ylim(0,24) +
  #stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = "Day length (hours") + 
  scale_color_manual(values = N_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position="none")
ggsave("output/Nouthern_dl.png", units = "cm" , height = 10, width = 10, dpi = 300)

################
# northern hemisphere without Fram
################
northern <- rbind(SPOT, L4, BBMO, BB) %>% 
  drop_na(Station) %>% 
  mutate(Week=as.integer(Week))

N_colors <- c("#56B4E9", "#006600", "#000066", "#E69F00")
ggplot(northern, aes(x=Week, y= dl, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 2) +
  ylim(7,17) +
  #stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = "Day length (hours") + 
  scale_color_manual(values = N_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position="none")
ggsave("output/Nouthern_no_FRAM_dl.png", units = "cm" , height = 10, width = 10, dpi = 300)
######################
## Back to southern
ggplot(southern, aes(x=Week, y= dl, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 2) +
  ylim(7,17) +
  #stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = "Day length (hours") + 
  scale_color_manual(values = S_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

ggsave("output/Southern_sst.png", units = "cm" , height = 10, width = 10, dpi = 300)

################
# Inorganic nitrogen
################
ggplot(northern, aes(x=Week, y= NOx, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 1.3) +
  stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = "NO2 + NO3 (umol L-1)") + 
  scale_color_manual(values = N_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position="none")
ggsave("output/Nouthern_NOx.png", units = "cm" , height = 10, width = 10, dpi = 300)


ggplot(southern, aes(x=Week, y= NOx, color=Station)) +
  theme_classic()+geom_point(shape=19, size = 1.5) +
  stat_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(y = "NO2 + NO3 (umol L-1)") + 
  scale_color_manual(values = S_colors) +
  theme(axis.title.x = element_text(size=14, vjust = 2),
        axis.title.y = element_text(size=14, vjust = 2),
        axis.text.y = element_text(size=14, vjust = 0.3, face="bold"),
        axis.text.x = element_text(size=14, vjust = 0.3, hjust=0.5, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position="none")

ggsave("output/Southern_NOx.png", units = "cm" , height = 10, width = 10, dpi = 300)
