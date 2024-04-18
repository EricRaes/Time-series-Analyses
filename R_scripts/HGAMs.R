library(mgcv)
library(gratia)
library(tidyverse)
library(tidygam)
library(forecast)
library(ggpubr)
library(patchwork)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading prefab files
# Create common dictionary for testing effects
# datetime, day_length, temp, chlorophyll, site (proxy for latitude also), year

# N.Hemisphere
FRAM <- read_csv("../Data_files/FRAM_Diversity_Meta.csv") %>% 
  mutate(site = "FRAM") %>% 
  mutate(Year = as.numeric(str_replace(Year, "Y_", ""))) %>%
  rename(latitude = lat) %>% 
  rename(longitude = lon) %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper, shannon) %>% 
  # Because of wildly different scales across sites, 
  #transforming predictors clarifies global relationship to diversity
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

Bedford <- read_csv("../Data_files/Bedford_Diversity_Meta.csv") %>% 
  mutate(site = "Bedford") %>% 
  rename(datetime = Date_Merge) %>% 
  rename(temp = Temperature) %>% 
  rename(chlorophyll = Chlorophyll.A) %>% 
  rename(longitude = longitude.x) %>% 
  rename(latitude = latitude.x)  %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper, shannon)  %>% 
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

L4 <- read_csv("../Data_files/L4_Diversity_Meta.csv")%>% 
  mutate(site = "L4") %>% 
  rename(day_length = day) %>% 
  rename(temp = Temperature..C.) %>% 
  rename(chlorophyll = Chlorophyll.A..ug.L.) %>% 
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper,shannon) %>% 
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

BBMO <- read_csv("../Data_files/BBMO_Diversity_Meta.csv")%>% 
  mutate(site = "BBMO") %>% 
  rename(day_length = ENV_Day_length_Hours_light) %>% 
  rename(temp = ENV_Temp) %>% 
  rename(chlorophyll = ENV_CHL_total) %>% 
  mutate(longitude = 2.48) %>% 
  mutate(latitude = 41.40) %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper,shannon) %>% 
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

SPOTS <- read_csv("../Data_files/SPOTS_Diversity_Meta.csv")%>% 
  mutate(site = "SPOTS") %>% 
  rename(day_length = daylength) %>% 
  rename(temp = sst) %>% 
  rename(chlorophyll = chla) %>% 
  mutate(longitude = 118.30) %>% 
  mutate(latitude = 33.30) %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper,shannon) %>% 
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

# S.Hemisphere
Maria <- read_csv("../Data_files/Maria_Diversity_Meta.csv")%>% 
  mutate(site = "Maria")%>% 
  #rename(day_length = daylength) %>% 
  rename(temp = CTDTemperature) %>% 
  rename(chlorophyll = CTDChlF_mgm3) %>%
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper,shannon) %>% 
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

Yongala <- read_csv("../Data_files/Yongala_Diversity_Meta.csv")%>% 
  mutate(site = "Yongala")%>% 
  #rename(day_length = daylength) %>% 
  rename(temp = CTDTemperature) %>% 
  rename(chlorophyll = CTDChlF_mgm3) %>% 
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper,shannon) %>% 
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

Rottnest <- read_csv("../Data_files/Rottnest_Diversity_Meta.csv")%>% 
  mutate(site = "Rottnest")%>% 
  #rename(day_length = daylength) %>% 
  rename(temp = CTDTemperature) %>% 
  rename(chlorophyll = CTDChlF_mgm3) %>% 
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  dplyr::select(datetime, day_length, temp, chlorophyll, 
                site, latitude, longitude,
                estimate, error, lower, upper,shannon) %>% 
  mutate(
    Daylength = (day_length), 
    Temperature = (temp), 
    Chlorophyll = (log(chlorophyll+1)))

# Bind
full <- Reduce(function(x, y) merge(x, y, all = TRUE), 
               list(FRAM, L4, Bedford, BBMO, SPOTS,
                    Maria, Yongala, Rottnest)) %>% 
  mutate(week = week(datetime)) %>% 
  mutate(month = month(datetime)) %>% 
  mutate(year = year(datetime))

full_na_omit <- na.omit(full)

pairs(full[,c("temp", 'chlorophyll', 'day_length', "week")])

############ Spell out model specifications

# Response variable = Diversity
# Covariates (fixed effects?) = Day length + Temperature + Chl-a
# Structure = Date?
# Random effects = site?

full <- full_na_omit %>% 
  arrange(site, datetime) %>%
  filter(Daylength < 20) %>% 
  filter(Daylength > 4) %>% 
  mutate(site = factor(site, levels = c("FRAM", "L4", "Bedford", "BBMO", "SPOTS",
                                        "Yongala", "Rottnest", "Maria"))) %>% 
  mutate()

################################################
# Daylength is not a measured variable like temp or chl. It therefore has infinite precision 
# and makes a less fair comparison for its relative effects on diversity. So, lets generate a new 
# column which is a random (normal) sample with the observed value as the mean and some sd < 1
# Test the GAMs additionally with Daylength_adj as a covariate

daylength_adj <- list()
DL <- full$Daylength

for(i in 1:length(DL)){
  random <- rnorm(1, mean=DL[[i]], sd=0.3)
  full$Daylength_adj[i] <- random

}

################################################
# Determining ARIMA structure for temporal autocorrelation
sites <- unique(full$site)

auto_sites <- list()

for (i in sites){
  auto_sites[[i]] <- auto.arima(full %>% 
                                filter(site == i) %>% 
                                  dplyr::select(estimate),
                                stepwise = F,
                                approximation = F)
  #paste(sites)
  #print(auto_sites)
}


#### GAM bulding #### 

hist((full$estimate), breaks = 20)
hist(log(full$estimate), breaks = 20)
hist((full$shannon), breaks = 20)


##########################################################
# Model S1 Parsimonious form, no environmental covs. 
##########################################################

mod_scale_pars <- mgcv::gam(log(estimate) ~ 
                            s(month, site, bs = "fs", k = 12)+
                            te(year, site, bs=c("tp", "re"))+
                            s(site, bs = "re"),
                          select = T,
                          #gamma = 1.4,
                          data=full, 
                          method="REML",
                          # KEY: Not Gaussian, which produces messy diagnostics
                          family=scat,
                          knots = list(month=c(0, 12)),
                          # Including autoregressive erros in model
                          correlation = corARMA(form = ~ 1|month, p = 3)
)

par(mfrow = c(2,2))
gam.check(mod_scale_pars)

png("../Figures/Model_pars_Diagnostic.png", width = 17, height = 12, units = "cm", res = 600)
par(mfrow = c(2,2))
appraise(mod_scale_pars)&theme_bw()
dev.off()

summary(mod_scale_pars)


par(mfrow = c(2,1))
plot(mod_scale_pars, shade = T)

# gratia approach
one <- draw(mod_scale_pars, residuals = F, select = 1)+
  geom_line(linewidth=1.15)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral",name = "Sites")+
  scale_x_continuous(limits = c(1,12), breaks = seq(1,12,1))+
  xlab("Month") +
  guides(color = guide_legend(override.aes = list(linewidth =3)))

two <- draw(mod_scale_pars, residuals = F, select = 2)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")

one + two + plot_layout(ncol = 1,
                      nrow = 2)
ggsave("../Figures/Model_parsimony_Partials.png", units = "cm", width = 17, height = 15, dpi = 600)

##########################################################
# Model G1 model form where global smoothers are specified, random intercepts are included for sites
##########################################################

# In case they were already run and saved (below), reload model here to avoid runtime
mod_scale_G <- read_rds("../Models/Model_G.rds")


mod_scale_G <- mgcv::gam(log(estimate) ~ 
                              s(Daylength_adj, k=15,bs = "tp",  m=2)+
                              s(Temperature, k=15,bs = "tp",  m=2)+
                              s(Chlorophyll, k=15,bs = "tp", m=2)+
                              #s(month, site , bs = "fs", k = 12)+
                              #te(year, month, bs=c("tp", "cc"))+
                           te(year, site, bs=c("tp", "re"))+
                           s(site, bs="re"),
                            select = T,
                            data=full, 
                            method="REML",
                            # KEY: Not Gaussian, which produces messy diagnostics
                            family=scat,
                            knots = list(month=c(1, 12)),
                            # Including autoregressive erros in model
                            correlation = corARMA(form = ~ 1|week, p = 3)
)



par(mfrow = c(2,2))
gam.check(mod_scale_G)
round(k.check(mod_scale_G),3)

png("../Figures/ModelG_Diagnostic.png", width = 17, height = 12, units = "cm", res = 600)
par(mfrow = c(2,2))
appraise(mod_scale_G)&theme_bw()
dev.off()

summary(mod_scale_G)

par(mfrow = c(2,2))
plot(mod_scale_G, shade = T)

# gratia approach
one <- draw(mod_scale_G, residuals = F, select = 1)+
  theme_bw(14)+
  ggtitle("")+
  scale_x_continuous(breaks = seq(8,18,1))+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")+
  xlab(expression(Temperatu))+
  xlab("Daylength (hours)")

two <- draw(mod_scale_G, residuals = F, select = 2)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")+
  xlab("Temperature (\u00B0C)")

three <- draw(mod_scale_G, residuals = F, select = 3)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")+
  xlab(expression(log(Chlorophyll~(mu*g~L^-1))))

four <- draw(mod_scale_G, residuals = F, select = 5)+
  theme_bw(14)+
  ggtitle("")+
  #scale_x_continuous(limits = c(1,12), breaks = seq(1,12,1))+
  scale_color_brewer(palette = "Spectral")


one + two + three + four + plot_layout(ncol = 2,
                                       nrow = 2)

ggsave("../Figures/Model_G_Partials.png", units = "cm", width = 17, height = 15, dpi = 600)



##########################################################
# S2 model form where a global smoother is not specified, but instead group level (site) is allowed
##########################################################

# In case they were already run and saved (below), reload model here to avoid runtime
mod_scale_S <- read_rds("../Models/Model_S2.rds")


mod_scale_S <- mgcv::gam(log(estimate) ~ 
                           s(Daylength_adj, site, bs = "fs",  m=2)+
                           s(Temperature, site , bs = "fs", m=2)+
                           s(Chlorophyll, site, bs = "fs", m=2)+
                           #s(month, site, bs="fs", m = 2)+
                           te(year, site, bs=c("tp", "re"))+
                           s(site, bs="re"),
                       select = T,
                       #gamma = 1.4,
                       data=full,
                       method="REML",
                       # KEY: Not Gaussian, which produces messy diagnostics
                       family=scat,
                       knots = list(month=c(0, 12)),
                       # Including autoregressive erros in model
                       correlation = corARMA(form = ~ 1|month, p = 3)
)

par(mfrow = c(2,2))
gam.check(mod_scale_S)
round(k.check(mod_scale_S),3)

#summary(mod_log)
summary(mod_scale_S)

png("../Figures/ModelS_Diagnostic.png", width = 17, height = 12, units = "cm", res = 600)
par(mfrow = c(2,2))
appraise(mod_scale_S)&theme_bw()
dev.off()

par(mfrow = c(2,3))
plot(mod_scale_S, shade = T)

AIC(mod_scale_G, mod_scale_pars, mod_scale_S)
BIC(mod_scale_G, mod_scale_pars, mod_scale_S)

# gratia approach
one <- draw(mod_scale_S, residuals = T, select = 1)+
  geom_line(linewidth=1.15)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  scale_x_continuous(breaks = seq(8,18,1))+
  xlab(expression(Temperatu))+
  xlab("Daylength (hours)")+
  theme(legend.position = "none")
  
two <- draw(mod_scale_S, residuals = T, select = 2)+
  geom_line(linewidth=1.15)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral", name="Sites")+
  xlab("Temperature (\u00B0C)")+
  guides(color = guide_legend(override.aes = list(linewidth =2.5)))
  

three <- draw(mod_scale_S, residuals = F, select = 3)+
  geom_line(linewidth=1.15)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  xlab(expression(log(Chlorophyll~(mu*g~L^-1))))+
  theme(legend.position = "none")
  

four <- draw(mod_scale_S, residuals = F, select = 5)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")
  scale_x_continuous(limits = c(1,12), breaks = seq(1,12,1))

five <- draw(mod_scale_S, residuals = F, select = 6)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")

one + two + three + four+ plot_layout(ncol = 2,
                                       nrow = 2)
ggsave("../Figures/Model_S_Partials.png", units = "cm", width = 17, height = 15, dpi = 600)

# Alternate approach, undesireable
draw(mod_scale_S, 
     select = c(1,2,3,5),
     residuals = T)&
  theme_bw(16)&
  ggtitle("")&
  scale_color_brewer(palette = "Spectral")&
  theme(legend.position = "none")


legend <- draw(mod_scale_S, residuals = T, select = 1)+
  theme_bw(18)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral", name = "Sites")+
  guides(color = guide_legend(override.aes = list(linewidth=2)))

legend <- as_ggplot(get_legend(legend))
ggsave("../Figures/Legend_Model_GS_Partials.png", units = "cm", width = 17, height = 15, dpi = 600)

##########################################################
# GS model form where a global smoother is specified by group level (site) is allowed
# If a site-level pattern deviates too far from global average it is penalized
# Computationally slow, minimal extra deviance explained, not useful
##########################################################

mod_scale_GS <- mgcv::gam(log(estimate) ~ 
                           s(Daylength, k=20, bs = "tp", m=2)+
                           s(Temperature, k=20, bs = "tp", m=2)+
                           s(Chlorophyll, k=20, bs = "tp",m=2)+
                           s(Daylength, site, bs = "fs",  m=1)+
                           s(Temperature, site , bs = "fs", m=1)+
                           s(Chlorophyll, site, bs = "fs", m=1)+
                            #s(month, site, bs="fs", m = 1)+
                            te(year, site, bs=c("tp", "re"), k=c(10, 10))+
                            s(site, bs = "re"),
                          select = T,
                          #gamma = 1.4,
                          data=full,
                          method="REML",
                          # KEY: Not Gaussian, which produces messy diagnostics
                          family=scat,
                          knots = list(month=c(0, 12)),
                          # Including autoregressive erros in model
                          correlation = corARMA(form = ~ 1|month, p = 3)
)

par(mfrow = c(2,2))
gam.check(mod_scale_GS)
round(k.check(mod_scale_GS),3)

summary(mod_scale_GS)

AIC(mod_scale_G, mod_scale_pars, mod_scale_GS, mod_scale_GS)
BIC(mod_scale_G, mod_scale_pars, mod_scale_GS, mod_scale_GS)


png("../Figures/ModelGS_Diagnostic.png", width = 17, height = 12, units = "cm", res = 600)
par(mfrow = c(2,2))
appraise(mod_scale_GS)&theme_bw()
dev.off()

par(mfrow = c(2,4))
plot(mod_scale_GS, shade = T)

# gratia approach
one <- draw(mod_scale_GS, residuals = F, select = 1)+
  geom_line(linewidth=1.1)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  scale_x_continuous(breaks = seq(8,18,1))+
  xlab(expression(Temperatu))+
  xlab("Daylength (hours)")+
  theme(legend.position = "none")

two <- draw(mod_scale_GS, residuals = F, select = 2)+
  geom_line(linewidth=1.1)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral", name="Sites")+
  xlab("Temperature (\u00B0C)")


three <- draw(mod_scale_GS, residuals = F, select = 3)+
  geom_line(linewidth=1.1)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  xlab(expression(log(Chlorophyll~(mu*g~L^-1))))+
  theme(legend.position = "none")


four <- draw(mod_scale_GS, residuals = F, select = 4)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")
#+
 # scale_x_continuous(limits = c(1,12), breaks = seq(1,12,1))

five <- draw(mod_scale_GS, residuals = F, select = 5)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")

six <- draw(mod_scale_GS, residuals = F, select = 6)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")

seven <- draw(mod_scale_GS, residuals = F, select = 8)+
  theme_bw(14)+
  ggtitle("")+
  scale_color_brewer(palette = "Spectral")+
  theme(legend.position = "none")

one + two + three + four + five + six + seven + plot_layout(ncol = 4,
                                        nrow = 2)

ggsave("../Figures/Model_GS_Partials.png", units = "cm", width = 17, height = 15, dpi = 600)


##################
# Save models
###################

write_rds(mod_scale_pars, file = "../Models/Model_S1_parsimony.rds")
write_rds(mod_scale_G, file = "../Models/Model_G.rds")
write_rds(mod_scale_S, file = "../Models/Model_S2.rds")
write_rds(mod_scale_GS, file = "../Models/Model_GS.rds")


