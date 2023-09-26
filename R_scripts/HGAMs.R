library(mgcv)
library(gratia)
library(tidyverse)
library(sjPlot)
library(gridExtra)
library(ggeffects)
library(tidygam)
library(forecast)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading prefab files
# Create common dictionary for testing effects
# datetime, day_length, temp, chlorophyll, site (proxy for latitude also), year

# Before standardizing predictors, I have log-transformed chlorophyll with a pseudo-count

# N.Hemisphere
FRAM <- read_csv("../Data_files/FRAM_Richness_Meta.csv") %>% 
  mutate(site = "FRAM") %>% 
  mutate(Year = as.numeric(str_replace(Year, "Y_", ""))) %>%
  rename(latitude = lat) %>% 
  rename(longitude = lon) %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper) %>% 
  # Because of wildly different scales across sites, 
  #transforming predictors clarifies global relationship to diversity
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

Bedford <- read_csv("../Data_files/Bedford_Richness_Meta.csv") %>% 
  mutate(site = "Bedford") %>% 
  rename(datetime = Date_Merge) %>% 
  rename(temp = Temperature) %>% 
  rename(chlorophyll = Chlorophyll.A) %>% 
  rename(longitude = longitude.x) %>% 
  rename(latitude = latitude.x)  %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper)  %>% 
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

L4 <- read_csv("../Data_files/L4_Richness_Meta.csv")%>% 
  mutate(site = "L4") %>% 
  rename(day_length = day) %>% 
  rename(temp = Temperature..C.) %>% 
  rename(chlorophyll = Chlorophyll.A..ug.L.) %>% 
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper) %>% 
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

BBMO <- read_csv("../Data_files/BBMO_Richness_Meta.csv")%>% 
  mutate(site = "BBMO") %>% 
  rename(day_length = ENV_Day_length_Hours_light) %>% 
  rename(temp = ENV_Temp) %>% 
  rename(chlorophyll = ENV_CHL_total) %>% 
  mutate(longitude = 2.48) %>% 
  mutate(latitude = 41.40) %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper) %>% 
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

SPOTS <- read_csv("../Data_files/SPOTS_Richness_Meta.csv")%>% 
  mutate(site = "SPOTS") %>% 
  rename(day_length = daylength) %>% 
  rename(temp = sst) %>% 
  rename(chlorophyll = chla) %>% 
  mutate(longitude = 118.30) %>% 
  mutate(latitude = 33.30) %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper) %>% 
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

# S.Hemisphere
Maria <- read_csv("../Data_files/Maria_Richness_Meta.csv")%>% 
  mutate(site = "Maria")%>% 
  #rename(day_length = daylength) %>% 
  rename(temp = CTDTemperature) %>% 
  rename(chlorophyll = CTDChlF_mgm3) %>%
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper) %>% 
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

Yongala <- read_csv("../Data_files/Yongala_Richness_Meta.csv")%>% 
  mutate(site = "Yongala")%>% 
  #rename(day_length = daylength) %>% 
  rename(temp = CTDTemperature) %>% 
  rename(chlorophyll = CTDChlF_mgm3) %>% 
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper) %>% 
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

Rottnest <- read_csv("../Data_files/Rottnest_Richness_Meta.csv")%>% 
  mutate(site = "Rottnest")%>% 
  #rename(day_length = daylength) %>% 
  rename(temp = CTDTemperature) %>% 
  rename(chlorophyll = CTDChlF_mgm3) %>% 
  rename(longitude = Longitude) %>% 
  rename(latitude = Latitude)  %>% 
  select(datetime, day_length, temp, chlorophyll, 
         site, latitude, longitude,
         estimate, error, lower, upper) %>% 
  mutate(
    estimate_scale_log = scale(log(estimate)),
    day_length = scale(day_length), 
         temp = scale(temp), 
         chlorophyll = scale(log(chlorophyll+1)))

# Bind
full <- Reduce(function(x, y) merge(x, y, all = TRUE), 
               list(FRAM, L4, Bedford, BBMO, SPOTS,
                    Maria, Yongala, Rottnest)) %>% 
  mutate(week = week(datetime)) %>% 
  mutate(month = month(datetime)) %>% 
  mutate(year = year(datetime))

full$site <- factor(full$site)
# Ordering causes issues with week random effect (?)
#full$site <- factor(full$site, levels = c("FRAM", "L4", "Bedford", "BBMO", "SPOTS",
#                                          "Maria", "Yongala", "Rottnest"))

full_na_omit <- na.omit(full)

pairs(full[,c("temp", 'chlorophyll', 'day_length', "week")])

############ Spell out model specifications

# Response variable = Richness
# Covariates (fixed effects?) = Day length + Temperature + Chl-a
# Structure = Date?
# Random effects = site?

full <- full %>% 
  arrange(site, datetime)

sites <- unique(full$site)

auto_sites <- list()

for (i in sites){
  auto_sites[[i]] <- auto.arima(full %>% 
                                filter(site == i) %>% 
                                select(estimate))
  #paste(sites)
  print(auto_sites)
}


# Test autocorrelation optima per time series

#### GAM bulding #### 

hist(log(full$estimate), breaks = 20)
#Nicer normality but significantly less variation explained
hist((full$estimate_scale_log), breaks = 20)

mod_log <- mgcv::gam(log(estimate) ~ s(day_length, bs = "cr", k = 15)+
             s(temp,bs = "tp", k = 15)+
             s(chlorophyll,bs = "tp", k = 15)+
             #s(month,bs = "cc", k = 15)+
             s(week,bs = "cc", k = 52)+
             s(site,bs = "re")+
             #s(year,bs = "tp", k = 5)+
             te(year,  site, bs=c("tp", "re"), k=c(10, 10), m=1) , 
           select = T,
           #gamma = 1.4,
           data=full, 
           method="REML",
           # KEY: Not Gaussian, which produces messy diagnostics
           family=scat,
           knots = list(week=c(0, 52)),
           # Including autoregressive erros in model
           correlation = corARMA(form = ~ 2|week, p = 1)
)

mod_scale <- mgcv::gam(estimate_scale_log ~ s(day_length, bs = "tp", k = 15)+
                       s(temp,bs = "tp", k = 15)+
                       s(chlorophyll,bs = "tp", k = 15)+
                         #te(day_length,  temp, bs=c("tp", "tp"), k=c(15, 15), m=1)+ 
                       #s(month,bs = "cc", k = 15)+
                       s(week,  bs = "cc", k = 52)+
                       s(site,bs = "re")+
                       #s(year,bs = "tp", k = 5)+
                       te(year,  site, bs=c("tp", "re"), k=c(10, 10), m=1) , 
                     select = T,
                     #gamma = 1.4,
                     data=full, 
                     method="REML",
                     # KEY: Not Gaussian, which produces messy diagnostics
                     family=scat,
                     knots = list(week=c(0, 52)),
                     # Including autoregressive erros in model
                     correlation = corARMA(form = ~ 1|week, q=3)
)

# High autocorrelation! Need to account
acf(resid(mod_log), main = "ACF", lag.max = 52)
pacf(resid(mod_scale), main = "pACF", lag.max = 52)

summary(mod_log)
summary(mod_scale)

# Log transform response (Diversity estimate) improves Q-Q normality and homoskedasticity
par(mfrow = c(2,2))
gam.check(mod_log)
gam.check(mod_scale)

par(mfrow = c(2,3))
plot(mod_log, shade = T)
plot(mod_scale, shade = T)

anova(mod_log, mod_scale)

plot(ggeffects::ggpredict(mod_auto), facets = T)
draw(mod_auto)

?plot_smooth
# Different marginal effecct approach
p <- plot_model(mod_auto, type = "pred", title  = "")

nCol <- floor(sqrt((length(p))))

do.call("grid.arrange", c(p, ncol=nCol))



