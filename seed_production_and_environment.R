# this file generates predicted fecundities based on plant trait measurements
# and analyzes relationships between fecundity and environment

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

# generate predicted fecundities
############
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

abundances <- read.csv("plec-and-valerianella-abundances-and-env-conditions.csv")

library(MASS)
# Best model for PLCO
P_m5 <- glm.nb(num_seeds ~ plant_height + number_of_inflors + infvol, 
               data = traits)
summary(P_m5)
P_m5.2 <- glm.nb(num_seeds ~ plant_height, 
                 data = traits)
summary(P_m5.2)

df_PLCO <- df_PLCO %>%
  mutate(num_seeds2 = 12.41195 + 1.03474*plant_height) 
## when using predict it was taking exp(coefficient*value) rather than exp(coefficient)*value
## so I manually computed seeds (num_seeds2) by multiply each value by exp(coefficient)

abundances$subplot_id<- 1:nrow(abundances) 
fitness_and_abundances_PLCO <- inner_join(df_PLCO, abundances)
fitness_and_abundances_PLCO[fitness_and_abundances_PLCO == "-"] <- NA
fitness_and_abundances_PLCO[fitness_and_abundances_PLCO == ""] <- NA

str(fitness_and_abundances_PLCO)


##########################
##########################

fitness_and_abundances <- read.csv("fitness_and_abundances.csv")
fitness_and_abundances_PLCO <- read.csv("fitness_and_abundances_PLCO.csv")
fitness_and_abundances_VALO <- read.csv("fitness_and_abundances_VALO.csv")


#################
# Grass Cover 

#m6 <- lme(num_seeds2 ~  as.numeric(X.grasscover.1m.2), 
#          data = fitness_and_abundances_PLCO, random=~1|Site,  na.action=na.omit) 
#summary(m6) # 
#m6.1 <- lme(num_seeds2 ~  as.numeric(X.grasscover.1m.2), 
#            data = fitness_and_abundances_VALO, random=~1|Site,  na.action=na.omit) 
#summary(m6.1) # 

m6.0 <- glmer.nb(num_seeds2 ~ as.numeric(X.grasscover.1m.2) + (1|Transect),
                 data = fitness_and_abundances_PLCO,  na.action=na.omit) 

m6 <- glmer.nb(num_seeds2 ~ as.numeric(X.grasscover.1m.2) + (1|Transect),
               data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m6.0))]),  na.action=na.omit) 
summary(m6) # (avg_soil_moisture_7.6cm) is not significant

m6.1 <- glmer.nb(num_seeds2 ~ (1|Transect), 
                 data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m6.0))]), na.action=na.omit) 
summary(m6.1)
anova(m6, m6.1, test = "LRT")

m7.0 <- glmer.nb(num_seeds2 ~ as.numeric(X.grasscover.1m.2) + (1|Transect),
                 data = fitness_and_abundances_VALO,  na.action=na.omit) 

m7 <- glmer.nb(num_seeds2 ~ as.numeric(X.grasscover.1m.2) + (1|Transect),
               data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m7.0))]),  na.action=na.omit) 
summary(m7) # (avg_soil_moisture_7.6cm) is not significant

m7.1 <- glmer.nb(num_seeds2 ~ (1|Transect), 
                 data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m7.0))]), na.action=na.omit) 
summary(m7.1)
anova(m7, m7.1, test = "LRT")


## newdat PLCO
newdat <- data.frame(X.grasscover.1m.2 <- seq(0, 100, 1))
newdat$y <- 35.79671 + 0.08754*as.numeric(newdat$X.grasscover.1m.2)  
newdat$Species = "Plectritis"
## newdat VALO
newdat2 <- data.frame(X.grasscover.1m.2 <- seq(0, 100, 1))
newdat2$y <- 32.83098 + 0.07710*as.numeric(newdat2$X.grasscover.1m.2)  
newdat2$Species = "Valerianella"

V <- ggplot(data = fitness_and_abundances, aes(x = as.numeric(X.grasscover.1m.2),
                                               y = num_seeds2))
V <- V + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
V <- V + geom_line(data = newdat, aes(x = X.grasscover.1m.2, y = y), size = 1)
V <- V + geom_line(data = newdat2, aes(x = X.grasscover.1m.2, y = y), size = 1, linetype= "dashed")
V <- V + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x=expression(paste("Percent grass cover (1 ", m^2, ")")), y="Predicted seeds per plant") +
  theme(legend.position="none")
V <- V + scale_color_manual(values=c("black", "black"))
V <- V + scale_shape_manual(values=c(19, 1))
V <- V + theme(axis.line = element_line(size = 2))
V <- V + theme(axis.title.x = element_text(vjust = 0,
                                           size = 20),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 20))
V <- V + theme(plot.title = element_text(size = 20, face = "bold"))
V <- V + theme(axis.line = element_line(colour = 'black', size = 1))
V




#######
# edit data frame for soil moisture. 

fitness_and_abundances_PLCO[fitness_and_abundances_PLCO == "#DIV/0!"] <- NA
fitness_and_abundances_VALO[fitness_and_abundances_VALO == "#DIV/0!"] <- NA
fitness_and_abundances[fitness_and_abundances == "#DIV/0!"] <- NA
#fitness_and_abundances_PLCO3 <- fitness_and_abundances_PLCO[!is.na(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm),]
fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm)
class(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm) 
fitness_and_abundances_VALO$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm)
class(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm) 
fitness_and_abundances$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances$avg_soil_moisture_7.6cm)
class(fitness_and_abundances$avg_soil_moisture_7.6cm) 



m8.0 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect),
               data = fitness_and_abundances_PLCO,  na.action=na.omit) 

m8 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect),
          data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m8.0))]),  na.action=na.omit) 
summary(m8) # (avg_soil_moisture_7.6cm) is not significant

m8.1 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m8.0))]), na.action=na.omit) 
summary(m8.1)
anova(m8, m8.1, test = "LRT")

m8.0.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect), 
                 data = fitness_and_abundances_VALO,  na.action=na.omit) 

m8.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_7.6cm + (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m8.0.2))]),  na.action=na.omit) 
summary(m8.2) # (soil.moisture.1.7.6cm) is not significant

m8.3 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m8.0.2))]),  na.action=na.omit) 
summary(m8.3) # (soil.moisture.1.7.6cm) is not significant
anova(m8.2, m8.3, test = "LRT")


X <- ggplot(data = fitness_and_abundances, aes(x = as.numeric(avg_soil_moisture_7.6cm),
                                               y = num_seeds2))
X <- X + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
X <- X + geom_segment(aes(x=0,xend=30,y=38,yend=38), color ="black", size = 1)
X <- X + geom_segment(aes(x=0,xend=30,y=34.59,yend=34.59), color ="black", size = 1, linetype = "dashed")
X <- X + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="Soil moisture 7.6 cm depth (%VWC)", y="Seeds / plant") +
  theme(legend.position="none")
X <- X + scale_color_manual(values=c("black", "black"))
X <- X + scale_shape_manual(values=c(19,1))
X <- X + theme(axis.line = element_line(size = 2))
X <- X + theme(axis.title.x = element_text(vjust = 0,
                                           size = 20),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 20))
X <- X + theme(plot.title = element_text(size = 20, face = "bold"))
X <- X + theme(axis.line = element_line(colour = 'black', size = 1))
X <- X + ggtitle("(A)")
X

# edit data frame further for deep soil moisture. 
fitness_and_abundances_PLCO$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances_PLCO$avg_soil_moisture_12cm)
class(fitness_and_abundances_PLCO$avg_soil_moisture_12cm) 
fitness_and_abundances_VALO$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances_VALO$avg_soil_moisture_12cm)
class(fitness_and_abundances_VALO$avg_soil_moisture_12cm) 
fitness_and_abundances$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances$avg_soil_moisture_12cm)
class(fitness_and_abundances$avg_soil_moisture_12cm) 

m9.0 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
          data = fitness_and_abundances_PLCO,  na.action=na.omit) 
summary(m9) # (avg_soil_moisture_7.6cm) is not significant

m9 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
                 data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m9.0))]),  na.action=na.omit) 
summary(m9) # (avg_soil_moisture_7.6cm) is not significant

m9.1 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_PLCO[ , all.vars(formula(m9.0))]),  na.action=na.omit) 
summary(m9.1)

anova(m9, m9.1, test = "LRT")

m9.0.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
                 data = fitness_and_abundances_VALO,  na.action=na.omit) 

m9.2 <- glmer.nb(num_seeds2 ~ avg_soil_moisture_12cm + (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m9.0.2))]),  na.action=na.omit) 
summary(m9.2) # (soil.moisture.12cm) is not significant


m9.3 <- glmer.nb(num_seeds2 ~ (1|Transect), 
            data = na.omit(fitness_and_abundances_VALO[ , all.vars(formula(m9.0.2))]),  na.action=na.omit) 
summary(m9.3) 

anova(m9.2, m9.3, test = "LRT")


Y <- ggplot(data = fitness_and_abundances, aes(x = as.numeric(avg_soil_moisture_12cm),
                                               y = num_seeds2))
Y <- Y + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
Y <- Y + geom_segment(aes(x=0,xend=40,y=38,yend=38), color ="black", size = 1)
Y <- Y + geom_segment(aes(x=0,xend=40,y=34.59,yend=34.59), color ="black", size = 1, linetype = "dashed")
Y <- Y + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="Soil moisture 12 cm depth (%VWC)", y="Seeds / plant") +
  theme(legend.position="none")
Y <- Y + scale_color_manual(values=c("black", "black"))
Y <- Y + scale_shape_manual(values=c(1, 19))
Y <- Y + theme(axis.line = element_line(size = 2))
Y <- Y + theme(axis.title.x = element_text(vjust = 0,
                                           size = 20),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 20))
Y <- Y + theme(plot.title = element_text(size = 20, face = "bold"))
Y <- Y + theme(axis.line = element_line(colour = 'black', size = 1))
Y <- Y + ggtitle("(B)")
Y


############
#soil depth

# soil with depth  less than 12 v greater than 12
fitness_and_abundances_temp2 <- fitness_and_abundances
PLCO_temp <- fitness_and_abundances_PLCO
VALO_temp <- fitness_and_abundances_VALO
fitness_and_abundances_temp2$soil.depth.1 <- revalue(fitness_and_abundances_temp2$soil.depth.1, c("<7"="<12cm", "7<x<12"="<12cm", ">12"=">12cm"))
PLCO_temp$soil.depth.1 <- revalue(PLCO_temp$soil.depth.1, c("<7"="<12cm", "7<x<12"="<12cm", ">12"=">12cm"))
VALO_temp$soil.depth.1 <- revalue(VALO_temp$soil.depth.1, c("<7"="<12cm", "7<x<12"="<12cm", ">12"=">12cm"))


m11.0 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
             data = PLCO_temp,  na.action=na.omit) 


m11.1 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
                  data = na.omit(PLCO_temp[ , all.vars(formula(m11.0))]),  na.action=na.omit) 
summary(m11.1)

m11.2 <- glmer.nb(num_seeds2 ~  (1|Transect), 
                  data = na.omit(PLCO_temp[ , all.vars(formula(m11.0))]),  na.action=na.omit) 
summary(m11.2)

anova(m11.1, m11.2, test = "LRT")

m11.0.2 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
                  data = VALO_temp,  na.action=na.omit) 

m11.3 <- glmer.nb(num_seeds2 ~ soil.depth.1 + (1|Transect), 
                  data = na.omit(VALO_temp[ , all.vars(formula(m11.0.2))]),  na.action=na.omit) 
summary(m11.3)

m11.4 <- glmer.nb(num_seeds2 ~ (1|Transect), 
                  data = na.omit(VALO_temp[ , all.vars(formula(m11.0.2))]),  na.action=na.omit) 
summary(m11.4)

anova(m11.3, m11.4, test = "LRT")


#order <- c("<7", "7<x<12", ">12")

Z <- ggplot(data=subset(fitness_and_abundances_temp2, !is.na(soil.depth.1)), aes(x=soil.depth.1, y=num_seeds2, color=Species))
Z <- Z + geom_boxplot(lwd=1)
#Z <- Z + 
#  scale_x_discrete(limits=order)
Z <- Z + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="soil depth", y="Seeds / plant") +
  theme(legend.position="none")
Z <- Z + scale_color_manual(values=c("black", "grey"))
Z <- Z + theme(axis.line = element_line(size = 2))
Z <- Z + theme(axis.title.x = element_text(vjust = 0,
                                           size = 20),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 20))
Z <- Z + theme(plot.title = element_text(size = 20, face = "bold"))
Z <- Z + theme(axis.line = element_line(colour = 'black', size = 1))
Z <- Z + ggtitle("(A)")
Z

############
############


require(gridExtra)
grid.arrange(X, Y, Z, ncol=1, nrow=3, 
             top=textGrob("", gp=gpar(fontsize=20,font=7))
)

ggarrange(#V + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
          Z + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
          X + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
          Y + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
          labels = c("(A)", "(B)", "(C)", "(D)"), font.label = list(size = 20))

