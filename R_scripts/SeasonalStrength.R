library(tidyverse)
library(devtools)
library(qiime2R)
library(ggplot2)
library(microbiome); packageVersion("microbiome")
library(magrittr)
library(breakaway)
library(phyloseq)
library(ggpubr)
library(timetk)
#library(ggfortify)
library(fpp3)
library(forecast)
library(seasonal)
library(tsfeatures)

##import data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

##########
########## clean up
########## 

carbom_no_mito = subset_taxa(carbom, !Order=="Chloroplast" )
carbom_no_mito = subset_taxa(carbom_no_mito, !Family=="Mitochondria" )
carbom_no_mito = subset_samples(carbom_no_mito, !Month== "XXX")
carbom_no_mito = subset_samples(carbom_no_mito, !Depth..m.== "60m")
carbom_no_mito = prune_samples(sample_sums(carbom_no_mito) >= 5000, carbom_no_mito)
carbom_no_mito = subset_samples(carbom_no_mito, sample_names(carbom_no_mito) != "BB14-22C")

#carbom_no_mito=prune_taxa(taxa_sums(carbom_no_mito)>=100, carbom_no_mito)

#carbom_no_mito=prune_samples(sample_data(carbom_no_mito)$Year=="2014", carbom_no_mito)

################################  
######################################## 
######################################## 

est_rich <- breakaway(carbom_no_mito, cutoff = 15)

#Some samples have extreme uncertainties based on the frequency counts
plot(est_rich, carbom_no_mito, color = "Year")

break_summary <- est_rich %>% summary
names(break_summary)[5] <- "SampleID"
META<- read.csv("METADATA_niskin_CTD_R_3.csv", header=TRUE)

brokenaway<-merge(META, break_summary, by="SampleID")
brokenaway$Depth..m. = factor(brokenaway$Depth..m., levels=c("1m", "5m", "10m"),
                              labels=c("euphotic","euphotic", "euphotic"))

brokenaway <- as_tibble(brokenaway)
brokenaway <- brokenaway[,c( 'SampleID', 'estimate', 'error', 'lower', 'upper',
                             'Year', 'Year.', 'month', 'day', 'Week','Month','ISO_Time', 'Depth..m.')]
brokenaway$datetime <-  as.Date(with(brokenaway, paste(Year, month, day,sep="-")), "%Y-%m-%d")
brokenaway <- brokenaway %>% arrange(datetime)

########
##plot it
p <- ggplot(data=brokenaway,aes(x=datetime, y=estimate, colour=Month))
p <- p + theme_bw(24)+
  geom_point(shape = 1,size = 2)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, size=1)+
  #  stat_regline_equation(aes(label = ..rr.label..), size=8)+
  geom_jitter(width = 0, height = 0)+
  scale_color_manual(values=c("blue1", "orange","#06f5fb",
                              "black" ,"red", "#a7cc7b",
                              "#cba9e5", "violet", "#808a8a",
                              "orange", "#ff249d", "black", "yellow" ))+
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

# Remove duplicates and group values by week (if plotted it fits original pattern quite well)
# The package (forecast) isn't great for weekly data
brokenaway_group <- brokenaway %>% 
  group_by(datetime) %>% 
  summarise_by_time(
    datetime, .by="month",
    estimate = median(estimate),
    error = median(error),
    lower = median(lower),
    upper = median(upper) 
)

# Index ts by week [1M] into tsibble object (packages cannot handle 1W data cleanly error prone)
brokenaway_group <- brokenaway_group %>%
  mutate(YearMonth = yearmonth(datetime)) %>%
  as_tsibble(index = YearMonth) %>% 
  # Impute NAs to missing values
  tsibble::fill_gaps()  

# where are the NAs?
which(is.na(brokenaway_group))
## Seasonal plots

brokenaway_group %>%
  gg_season(estimate, labels = "both") +
  labs(y = "Richness",
       title = "Bedford Basin")

# Seasonal subseries just to play around
brokenaway_group %>%
  gg_subseries(estimate) +
  labs(y = "Richness",
       title = "Bedford Basin")
# Blue h-line indicates monthly average

# lag plot 
brokenaway_group %>%
  gg_lag(estimate, geom = "point") +
  labs(y = "Richness",
       title = "Bedford Basin")
# strong lag, especially lag1-2 indicating pronounced seasonality

#autocorrelation
brokenaway_group %>%
  # Max for # of weeks
  ACF(estimate, lag_max = 52) %>%
  autoplot() + 
  labs(
       title = "Bedford Basin")
# unsurprisingly, data closest together (in 1W increments) are highly correlated and become much more similar as winter sets in again

#x11 method
x11_dcmp <- brokenaway_group %>%
  model(x11 = X_13ARIMA_SEATS(estimate ~ x11())) %>%
  components()

autoplot(x11_dcmp) +
  labs(title =
         "Decomposition of sequence richness using X-11.")

x11_dcmp %>%
  ggplot(aes(x = YearMonth)) +
  geom_line(aes(y = estimate, colour = "Data")) +
  geom_line(aes(y = season_adjust,
                colour = "Seasonally Adjusted")) +
  geom_line(aes(y = trend, colour = "Trend")) +
  labs(y = "Richness",
       title = "Bedford Basin") +
  scale_colour_manual(
    values = c("gray", "#0072B2", "#D55E00"),
    breaks = c("Data", "Seasonally Adjusted", "Trend")
  )


# SEATS method...
seats_dcmp <- brokenaway_group %>%
  model(seats = X_13ARIMA_SEATS(estimate ~ seats())) %>%
  components()
autoplot(seats_dcmp) +
  labs(title =
         "Decomposition of sequence richness using SEATS")
# not very good at pulling seasonality

# STL is more robust in other fields and based on my limited knowledge, I would recommend it
brokenaway_group %>%
  model(
    STL(estimate ~ trend(window = 13) +
          season(window = "periodic"),
        robust = TRUE)) %>%
  components() %>%
  autoplot()
# This is a bit better, need to play with trend and season window 

brokenaway_group %>%
  features(estimate, feat_stl) %>% 
  ggplot(aes(x = trend_strength, y = seasonal_strength_year)) +
  geom_point(size=5)
# the seaonal strength is much greater than trend, add below argument (plus color) when other time series are included
  facet_wrap((station))
 
# Ugly plot it has one point (Bedford Basin), but if you include the five other timeseries, you can directly 
# compare which have the strongest seasonal component and whether this is latitudinally structured

  # Here is the punchline: the variance in seasonality should outweigh trend and random variation,
  # The relative strength of seasonaility can be quantified and comapred across time series!
  
# shortcut way to pull ts features
myfeatures <- tsfeatures(brokenaway_group)



