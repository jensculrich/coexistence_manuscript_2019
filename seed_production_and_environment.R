# this file contains code to generate predicted 
# fecundities based on plant trait measurements
# and to analyze relationships between fecundity 
# and soil moisture and soil depth

library(tidyverse)
library(stats4)
library(nlme)
library(gridExtra)
library(grid)
library(broom)
library(purrr)
library(repurrrsive)
library(listviewer)
library(varhandle)
library(plyr)
library(lme4)
library(MASS)

###############################################################################
# select best model to predict seed production of field plants                #
# based on allometric measurements and seed counts on a sample of individuals #
###############################################################################

df <- read.csv("fitness_proxies_from_the_field_2.csv") #read data file

df2 <- df %>% # reshape data into long form with key representing type of measurement and value the measurement itself
  gather(key, value, -Site, -Area, -Transect, -Plot, -Subplot, -Date, -Species)

### add a number (called plant_id) 1:2760 x 4.
plant_id <- rep(1:2760, times=4, each=1)
df3 <- as.data.frame(cbind(df2, plant_id))
### spread by key to have a list of measurements taken from each of 2760 plants
# rename key, all ht1, ht2, ht3 etc need to be relabelled as plain ht
# same for all other key types
# just make a new column with this sequence: 2670 as ht, 2670 as x.inf, 2670 as inf_ht, 2670 as inf_w
new_key <- rep(c("ht", "x.inf", "inf_ht", "inf_w"), times=1, each=2760)
df3 <- as.data.frame(cbind(df3, new_key)) #cbind with data frame
df4 <- df3[ , -8] # drop old key  

df5 <- df4 %>% # reshape data into long form with key representing type of measurement and value the measurement itself
  spread(key=new_key, value=value)

# add a subplot ID to identify unique subplots (there are 20 rows per subplot)
subplot_id <- rep(1:138, times=1, each=20)
df6 <- as.data.frame(cbind(df5, subplot_id)) #cbind with data frame

# rename columns to be consistent with model parameters
colnames(df6)[which(names(df6) == "ht")] <- "plant_height"
colnames(df6)[which(names(df6) == "x.inf")] <- "number_of_inflors"
colnames(df6)[which(names(df6) == "inf_ht")] <- "inflor_height"
colnames(df6)[which(names(df6) == "inf_w")] <- "inflor_width"

str(df6) #plant trait measures are factors
#convert to numeric
df6[, 9:12] <- sapply(df6[, 9:12], as.numeric)
str(df6)

# calculate volume as cylinder if plectritis (pi*r*2h)
df_PLCO <- df6 %>% 
  filter(Species == "PLCO") %>% # filter by PLCO
  mutate(infvol=pi*inflor_width*inflor_height) # *1/2 width (to get radius) and *2 to get height cancel out

# and volume as hemisphere if valerianella ((2 / 3)*pi*r3)
df_VALO <- df6 %>% 
  filter(Species == "VALO") %>% # filter by VALO 
  mutate(infvol=(2/3)*pi*(1/2*inflor_width)^3)

df7 <- as.data.frame(rbind(df_PLCO, df_VALO))
head(df7)
str(df7)
View(df7)

# Now use new df7 (individual plant measurements kept to explore effects on fitness)
# First predict seed production
traits <- read.csv("cowichan_pollen_supplementations.csv") 
View(traits)

colnames(traits)[which(names(traits) == "height")] <- "plant_height"
colnames(traits)[which(names(traits) == "num_infl")] <- "number_of_inflors"
colnames(traits)[which(names(traits) == "inf_volume")] <- "infvol"
colnames(traits)[which(names(traits) == "Fruit")] <- "num_seeds"

# Best model for PLCO
P_m5 <- glm.nb(num_seeds ~ plant_height + number_of_inflors + infvol, 
               data = traits)
summary(P_m5)
P_m5.2 <- glm.nb(num_seeds ~ plant_height, 
                 data = traits)
summary(P_m5.2)

df_PLCO <- df_PLCO %>%
  mutate(num_seeds2 = 12.41195 + 1.03474*plant_height) 

abundances <- read.csv("plec-and-valerianella-abundances-and-env-conditions.csv")
abundances$subplot_id<- 1:nrow(abundances) 

fitness_and_abundances_PLCO <- inner_join(df_PLCO, abundances)
fitness_and_abundances_PLCO[fitness_and_abundances_PLCO == "-"] <- NA
fitness_and_abundances_PLCO[fitness_and_abundances_PLCO == ""] <- NA

str(fitness_and_abundances_PLCO)


#############################################
# Relationships between soil moisture/depth #
# and predicted seed production.            #
#############################################

fitness_and_abundances <- read.csv("fitness_and_abundances.csv")
fitness_and_abundances_PLCO <- read.csv("fitness_and_abundances_PLCO.csv")
fitness_and_abundances_VALO <- read.csv("fitness_and_abundances_VALO.csv")

#######
# edit data frame for soil moisture. 

fitness_and_abundances_PLCO[fitness_and_abundances_PLCO == "#DIV/0!"] <- NA
fitness_and_abundances_VALO[fitness_and_abundances_VALO == "#DIV/0!"] <- NA
fitness_and_abundances[fitness_and_abundances == "#DIV/0!"] <- NA
fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm)
class(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm) 
fitness_and_abundances_VALO$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm)
class(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm) 
fitness_and_abundances$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances$avg_soil_moisture_7.6cm)
class(fitness_and_abundances$avg_soil_moisture_7.6cm) 

# SEED PRODUCTION V SOIL MOISTURE (7.6 cm depth) with random effect of transect
# Plectritis Response
m8.0 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect),
               data = fitness_and_abundances_PLCO,  na.action=na.omit) 

m8 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect),
          data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m8.0))]),  na.action=na.omit) 
summary(m8)

m8.1 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m8.0))]), na.action=na.omit) 
summary(m8.1)
anova(m8, m8.1, test = "LRT")

# Valerianella Response
m8.0.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect), 
                 data = fitness_and_abundances_VALO,  na.action=na.omit) 

m8.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m8.0.2))]),  na.action=na.omit) 
summary(m8.2)

m8.3 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m8.0.2))]),  na.action=na.omit) 
summary(m8.3)
anova(m8.2, m8.3, test = "LRT")


# edit data frame further to use deep soil moisture as a fixed effect. 
fitness_and_abundances_PLCO$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances_PLCO$avg_soil_moisture_12cm)
class(fitness_and_abundances_PLCO$avg_soil_moisture_12cm) 
fitness_and_abundances_VALO$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances_VALO$avg_soil_moisture_12cm)
class(fitness_and_abundances_VALO$avg_soil_moisture_12cm) 
fitness_and_abundances$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances$avg_soil_moisture_12cm)
class(fitness_and_abundances$avg_soil_moisture_12cm) 

# SEED PRODUCTION V SOIL MOISTURE (12 cm depth) with random effect of transect
# Plectritis Response
m9.0 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
          data = fitness_and_abundances_PLCO,  na.action=na.omit) 
summary(m9) 

m9 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
                 data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m9.0))]),  na.action=na.omit) 
summary(m9)

m9.1 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m9.0))]),  na.action=na.omit) 
summary(m9.1)

anova(m9, m9.1, test = "LRT")

# Valerianella response
m9.0.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
                 data = fitness_and_abundances_VALO,  na.action=na.omit) 

m9.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m9.0.2))]),  na.action=na.omit) 
summary(m9.2) 

m9.3 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m9.0.2))]),  na.action=na.omit) 
summary(m9.3) 

anova(m9.2, m9.3, test = "LRT")

############
#soil depth

# edit data frame to classify soil with depth less than 12 v greater than 12 cm
fitness_and_abundances_temp2 <- fitness_and_abundances
PLCO_temp <- fitness_and_abundances_PLCO
VALO_temp <- fitness_and_abundances_VALO
fitness_and_abundances_temp2$soil.depth.1 <- revalue(fitness_and_abundances_temp2$soil.depth.1, c("<7"="<12cm", "7<x<12"="<12cm", ">12"=">12cm"))
PLCO_temp$soil.depth.1 <- revalue(PLCO_temp$soil.depth.1, c("<7"="<12cm", "7<x<12"="<12cm", ">12"=">12cm"))
VALO_temp$soil.depth.1 <- revalue(VALO_temp$soil.depth.1, c("<7"="<12cm", "7<x<12"="<12cm", ">12"=">12cm"))


# SEED PRODUCTION V SOIL DEPTH (12 cm depth) with random effect of transect
# Plectritis Response
m11.0 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
             data = PLCO_temp,  na.action=na.omit) 


m11.1 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
                  data = na.omit(PLCO_temp[ , all.vars(formula(m11.0))]),  na.action=na.omit) 
summary(m11.1)

m11.2 <- glmer.nb(num_seeds2 ~  (1|Transect), 
                  data = na.omit(PLCO_temp[ , all.vars(formula(m11.0))]),  na.action=na.omit) 
summary(m11.2)

anova(m11.1, m11.2, test = "LRT")

# Valerianella response
m11.0.2 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
                  data = VALO_temp,  na.action=na.omit) 

m11.3 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
                  data = na.omit(VALO_temp[ , all.vars(formula(m11.0.2))]),  na.action=na.omit) 
summary(m11.3)

m11.4 <- glmer.nb(num_seeds2 ~ (1|Transect), 
                  data = na.omit(VALO_temp[ , all.vars(formula(m11.0.2))]),  na.action=na.omit) 
summary(m11.4)

anova(m11.3, m11.4, test = "LRT")