library(tidyr )
library(vegan)
library(scales)
library(grid)
library(reshape2)
library(ggplot2)
library(QsRutils)
library(ggplot2)
library(data.table)
library(ggpubr)
library("plyr")
library(dplyr)
library(phyloseq)
library(knitr)
library(tibble)
library(ggrepel)
library(BiocManager)
library(microbiome)
library(aweek)
library(knitr)
library(ggords)
library(tidyverse)
library(ggrepel)
library(iNEXT)
library(pairwiseAdonis)
library(eulerr)


######## Import data as an example here we show the Bedford Basin data
otu_mat <-"ASV_Table_Bedford_Basin.csv"
tax_mat <- "Taxonomy_species_Silva138.1_Bedford_Basin.csv"
samples_df <- "METADATA_Bedford_Basin.csv"
file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

Bedford<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

######## Remove taxa
Bedford_no_mito = subset_taxa(Bedford, !Kingdom=="NA")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Kingdom=="Eukaryota")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Phylum=="NA")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Order=="Chloroplast")
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria")
Bedford_no_mito = subset_samples(Bedford_no_mito, Depth < 50) # remove deep samples

############ have a look at the read distribution
#Make a data frame with a column for the read counts of each sample##### 
sample_sum_df <- data.frame(sum = sample_sums(Bedford_no_mito))
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth - Bedford Basin 16S rRNA ASV data") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

########### rarefaction curves
library(devtools)
devtools::install_github("gauravsk/ranacapa")

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}
p <- ggrare(Bedford_no_mito, step = 10, color = "Month", se = FALSE)

########## Set factor levels
sample_data(Bedford_no_mito)$Month = factor(sample_data(Bedford_no_mito)$Month, levels=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sep","Oct", "Nov", "Dec"))
sample_data(Bedford_no_mito)$Year. = factor(sample_data(Bedford_no_mito)$Year., levels=c("Year_2014", "Year_2015", "Year_2016", "Year_2017"))

######### Set colours
my.cols <- c("darkorchid","orange", "blue1",
             "green" ,"black", 
             "cyan", "red", "saddlebrown",
             "goldenrod", "violet",
             "gold1", "slateblue", "chartreuse")

#########################
#Ordination PCA CLR and Euclidian distance
#######################
set.seed(66)
Bedford_clr <- microbiome::transform(Bedford_no_mito, "clr")   
out.pcoa.logt <- ordinate(Bedford_clr, method = "RDA", distance = "euclidean")
evals <- out.pcoa.logt$CA$eig
p3<-plot_ordination(Bedford_clr, out.pcoa.logt, type = "Sample", 
                    color = "Month" ) 
#p3$layers <- p3$layers[-1]
p3+ggtitle("PCA Bedford 16S rRNA ASV")+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.7))+
  theme(legend.key=element_blank())+
  geom_point(size=3)+
  scale_colour_manual(values=my.cols)

######################################
#Alpha diversity and Correlations     #change sample.size = 5000, to your read depth of interest
########################################

Bedford_rar <-rarefy_even_depth(Bedford_no_mito, sample.size = 5000,
                               rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)

results = estimate_richness(Bedford_rar, measures =c( 'Chao1', "Shannon"))
write.csv(results, "ASV_alpha.csv")
Richness<-read.csv("ASV_alpha.csv", header=TRUE)
names(Richness)[1] <- "X"
Richness$X <- gsub("[.]", "-", Richness$X)

d = sample_data(Bedford_rar)
d1<- data.frame(d)
write.csv(d1, "META_ASV_alpha.csv")
META<-read.csv("META_ASV_alpha.csv", header=TRUE)
#names(write.csv)[1] <- "X"

METADATA<-merge(Richness,META, by="X", all=TRUE)

#changing name of first column as dates
names(METADATA)[5] <- "Date_Merge"
#checking class of dates
class(METADATA$Date_Merge)
#converting as dataframe
METADATA$Date_Merge <- as.Date(METADATA$Date_Merge, format = "%d/%m/%Y") #%Y/%m/%d - %d/%m/%Y
#checking class again
class(METADATA$Date_Merge)

### stats for alpha
Chao1_Means<-compare_means(Chao1 ~ Month,  data = METADATA, p.adjust.method = "bonferroni")
Shannon_Means<-compare_means(Shannon ~ Month,  data = METADATA, p.adjust.method = "bonferroni")
ddply(METADATA, ~Year, plyr:::summarise, mean = mean(Chao1), sd = sd(Chao1),
      max = max(Chao1), min = min(Chao1))

### Plotting correlations

## Set factor levels
METADATA$Depth..m. = factor(METADATA$Depth..m., levels=c("1m", "5m", "10m"), 
                            labels=c("euphotic","euphotic", "euphotic"))
METADATA$Month = factor(METADATA$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

## Set colours
my.cols1 <- c("darkorchid","orange", "blue1",
              "green" ,"black", 
              "cyan", "red", "saddlebrown",
              "goldenrod", "violet",
              "gold1", "slateblue", "chartreuse")

## ggplot
p<-ggplot(data=METADATA, aes(x=day_length , y=Chao1,  fill=Depth..m.)) 
# use #day_length #Chlorophyll.A #Temperature
p <- p + geom_point(aes(color=Depth..m.))+
  geom_point(aes(color=Depth..m.),shape = 21,size = 2,colour = "black")+
  scale_fill_manual(values=my.cols1) +
  scale_colour_manual(values=my.cols1)+
  theme(axis.title.x = element_text(size=16, vjust = 0),
        axis.title.y = element_text(size=16, vjust = 0),
        axis.text.y = element_text(size=22, vjust = 0.5),
        axis.text.x = element_text(size=22, vjust = 0, angle = 0, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Day length")+ylab("Diversity")+ggtitle("...")+
  scale_x_continuous(limits=c(8, 16), breaks=seq(8,16,2))

p
p+stat_smooth(method=lm)+
  stat_cor(method = "pearson", aes(color = Depth..m.), size=6, label.x = 12, p.accuracy = 0.001, r.accuracy = 0.01)+
  guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)

###################################
# Alpha Diversity and time
################################

Bedford_rar <-rarefy_even_depth(Bedford_no_mito, sample.size = 5000,
                               rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)

### for the rare microbiome, ASVs contributing less < 1% 
## if you rarefy to  5000 reads, then look at ASVs which contribute less than 50 reads. 50 is 1% of 5000 reads
## if you rarefy to 80000 reads, then change 50 to 800
#Bedford_rar= prune_taxa(taxa_sums(Bedford_rar) < 50, Bedford_rar)

sample_data(Bedford_rar)$Depth..m. = factor(sample_data(Bedford_rar)$Depth..m.,levels=c("1m", "5m", "10m"), 
                                           labels=c("euphotic","euphotic", "euphotic"))

p<-plot_richness(Bedford_rar,x="Week", color="Year.", measures=c("Chao1", "Shannon")) #Depth..m.
p + theme_bw()+geom_point(shape = 1,size = 2,colour = "black")+
  scale_color_manual(values=c("blue1", "orange","#06f5fb",
                              "black" ,"red", "#a7cc7b",
                              "#cba9e5", "violet", "#808a8a",
                              "orange", "#ff249d" ))+
  ggtitle("Bedford Basin") +
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=14, vjust = 0.3, angle = 0, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"))+
  theme(panel.grid.major = element_blank())+
  stat_smooth(method = "loess",  aes(color = Depth..m.), formula = y ~ x, se = TRUE)

######################################























