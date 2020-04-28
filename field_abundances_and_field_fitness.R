# this file contains code to analyze relationships between fecundity 
# and hetero-/con-specific densities
# and between abundance and hetero/conspecific densities

library(tidyverse)
library(stats4)
library(nlme)
library(gridExtra)
library(grid)
library(broom)
library(purrr)
library(repurrrsive)
library(listviewer)
library(ggpubr)
library(MASS)
library(lme4)


fitness_and_abundances <- read.csv("fitness_and_abundances.csv")
fitness_and_abundances <- fitness_and_abundances %>%
  mutate(Transect = factor(Transect)) %>%
  mutate(Plot = factor(Plot))

fitness_and_abundances_PLCO <- read.csv("fitness_and_abundances_PLCO.csv")
fitness_and_abundances_PLCO <- fitness_and_abundances_PLCO %>%
  mutate(Transect = factor(Transect)) %>%
  mutate(Plot = factor(Plot))

fitness_and_abundances_VALO <- read.csv("fitness_and_abundances_VALO.csv")
fitness_and_abundances_VALO <- fitness_and_abundances_VALO %>%
  mutate(Transect = factor(Transect)) %>%
  mutate(Plot = factor(Plot))

### remove the big outlier from the data set (many more seeds preditced vs any other plot)
which.max(fitness_and_abundances_PLCO$num_seeds2)
fitness_and_abundances_PLCO2 <- fitness_and_abundances_PLCO[-31, ]

### remove the big outlier from the data set (more than twice as many Plectritis density than second highest value)
which.max(fitness_and_abundances_VALO$X.Plectritis..including.pollinator.focal.plants.)
fitness_and_abundances_VALO2 <- fitness_and_abundances_VALO[-9, ]

### remove thetis lake data, where no Valerianella found
fitness_and_abundances_VALO2_filtered <- fitness_and_abundances_VALO2 %>% 
  filter(Site != "Thetis")
fitness_and_abundances_PLCO2_filtered <- fitness_and_abundances_PLCO2 %>% 
  filter(Site != "Thetis")


############################################################
# Relationships between fecundity and competitor densities #
############################################################

# plectritis fitness against plectritis density with mixed effects glm.nb 
str(fitness_and_abundances_PLCO2_filtered)

mixed_m1 <- glmer.nb(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants. 
                     + (1|Transect/Plot), 
          data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m1)

## estimate intercept without independent variables
mixed_m2 <- glmer.nb(num_seeds2 ~ (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m2)
anova(mixed_m1, mixed_m2, test = "LRT")

# valerianella fitness against plectritis density with mixed effects glm.nb 
mixed_m3 <- glmer.nb(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants. 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_VALO2_filtered, na.action=na.omit) 
summary(mixed_m3)

mixed_m4 <- glmer.nb(num_seeds2 ~ 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_VALO2_filtered, na.action=na.omit) 
summary(mixed_m4)
anova(mixed_m4)
anova(mixed_m3, mixed_m4, test = "LRT")

#################
## plectritis fitness against valerianella density with mixed effects glm.nb 
mixed_m5 <- glmer.nb(num_seeds2 ~ X.Valerianella 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m5) 

# estimate intercept without independent variables
mixed_m6 <- glmer.nb(num_seeds2 ~ (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m6) 
anova(mixed_m5, mixed_m6, test = "LRT")

# plectritis fitness against valerianella density with mixed effects glm.nb 
mixed_m7 <- glmer.nb(num_seeds2 ~ X.Valerianella 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_VALO2_filtered,  na.action=na.omit) 
summary(mixed_m7) 

mixed_m8 <- glmer.nb(num_seeds2 ~ (1|Transect/Plot), 
                     data = fitness_and_abundances_VALO2_filtered,  na.action=na.omit) 
summary(mixed_m8) 
anova(mixed_m7, mixed_m8, test = "LRT")


###################################################
# analyze relationship between species abundances #
###################################################

df <- read.csv("plec-and-valerianella-abundances-and-env-conditions.csv")

View(df)
df[df == "-"] <- NA
df[df == ""] <- NA
df[df == "#DIV/0!"] <- NA


# first filter out Thetis lake data 
df_filtered <- df %>% 
  filter(Site != "Thetis")

# relationship between abundances at 0.1 m^2 density
mixed_m9 <- glmer.nb(X.Valerianella ~ X.Plectritis..including.pollinator.focal.plants. 
                     + (1|Site/Transect/Plot), 
                     data = df_filtered,  na.action=na.omit) 
summary(mixed_m9) 

mixed_m9.2 <- glmer.nb(X.Valerianella ~  
                     + (1|Site/Transect/Plot), 
                     data = df_filtered,  na.action=na.omit) 
summary(mixed_m9.2) 

anova(mixed_m9, mixed_m9.2, test = "LRT")

##########
##########

# relationship between abundances at 1m^2 density

mixed_m10 <- glmer.nb(X.Valerianella.1m.2..including.subplots. ~ X.Plectritis.1m.2..including.subplots. 
                     + (1|Site/Transect/Plot), 
                     data = df_filtered,  na.action=na.omit) 
summary(mixed_m10)

mixed_m10.2 <- glmer.nb(X.Valerianella.1m.2..including.subplots. ~  
                      + (1|Site/Transect/Plot), 
                      data = df_filtered,  na.action=na.omit) 
summary(mixed_m10.2)

anova(mixed_m10, mixed_m10.2, test = "LRT")
