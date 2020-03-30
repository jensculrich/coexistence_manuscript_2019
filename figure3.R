# figure3

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
library(MASS)

##########################
##########################

fitness_and_abundances <- read.csv("fitness_and_abundances.csv")
fitness_and_abundances_PLCO <- read.csv("fitness_and_abundances_PLCO.csv")
fitness_and_abundances_VALO <- read.csv("fitness_and_abundances_VALO.csv")

fitness_and_abundances_PLCO[fitness_and_abundances_PLCO == "#DIV/0!"] <- NA
fitness_and_abundances_VALO[fitness_and_abundances_VALO == "#DIV/0!"] <- NA
fitness_and_abundances[fitness_and_abundances == "#DIV/0!"] <- NA
fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm)
class(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm) 
fitness_and_abundances_VALO$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm)
class(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm) 
fitness_and_abundances$avg_soil_moisture_7.6cm <- unfactor(fitness_and_abundances$avg_soil_moisture_7.6cm)
class(fitness_and_abundances$avg_soil_moisture_7.6cm) 


##########
#mixed_m1 <- glmer.nb(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants. 
#                     + (1|Transect/Plot), 
#                     data = fitness_and_abundances_PLCO_filtered,  na.action=na.omit) 
#summary(mixed_m1) # singular boundary due to small capture ofvarianve by random effects
# Plectritis is non-significant
#anova(mixed_m1)

#################
# Grass Cover 

m6 <- glm.nb(num_seeds2 ~  as.numeric(X.grasscover.1m.2), 
          data = fitness_and_abundances_PLCO,  na.action=na.omit) 
summary(m6) 


m6.1 <- glm.nb(num_seeds2 ~  as.numeric(X.grasscover.1m.2), 
            data = fitness_and_abundances_VALO,  na.action=na.omit) 
summary(m6.1) # X.grasscover.1m.2 is significant (positive)
## intercept is significant = 34.82207
anova(m6)
anova(m6.1)

newdata <- data.frame(
  X.grasscover.1m.2 = 
    rep(seq(from = min(fitness_and_abundances_PLCO$X.grasscover.1m.2, na.rm = TRUE), 
            to = max(fitness_and_abundances_PLCO$X.grasscover.1m.2, na.rm = TRUE), length.out = 100)))
newdata <- cbind(newdata, predict(m6, newdata, type = "link", se.fit=TRUE))
newdata <- within(newdata, {
  num_seeds2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})



V <- ggplot(data = fitness_and_abundances, aes(x = as.numeric(X.grasscover.1m.2),
                                               y = num_seeds2))
V <- V + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
V <- V + geom_line(data = newdata, aes(x = X.grasscover.1m.2, 
                                       y = num_seeds2), size = 1, colour = "black")

V <- V + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(),panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x=expression(paste("Percent grass cover (1 ", m^2, ")")), y="Seeds / plant") +
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
V <- V + geom_ribbon(data = newdata, aes(ymin = LL, ymax = UL), alpha = .25)
V

m8 <- glm.nb(num_seeds2 ~ avg_soil_moisture_7.6cm,
               data = fitness_and_abundances_PLCO,  na.action=na.omit) 
summary(m8) # (avg_soil_moisture_7.6cm) is not significant
## intercept is significant = exp(3.417535)
anova(m8)

newdata2 <- data.frame(
  avg_soil_moisture_7.6cm = 
    rep(seq(from = min(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm, na.rm = TRUE), 
            to = max(fitness_and_abundances_PLCO$avg_soil_moisture_7.6cm, na.rm = TRUE), length.out = 100)))
newdata2 <- cbind(newdata2, predict(m8, newdata2, type = "link", se.fit=TRUE))
newdata2 <- within(newdata2, {
  num_seeds2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

m8.2 <- glm.nb(num_seeds2 ~ avg_soil_moisture_7.6cm, 
            data = fitness_and_abundances_VALO,  na.action=na.omit) 
summary(m8.2)

newdata3 <- data.frame(
  avg_soil_moisture_7.6cm = 
    rep(seq(from = min(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm, na.rm = TRUE), 
            to = max(fitness_and_abundances_VALO$avg_soil_moisture_7.6cm, na.rm = TRUE), length.out = 100)))
newdata3 <- cbind(newdata3, predict(m8.2, newdata2, type = "link", se.fit=TRUE))
newdata3 <- within(newdata3, {
  num_seeds2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


X <- ggplot(data = fitness_and_abundances, aes(x = as.numeric(avg_soil_moisture_7.6cm),
                                               y = num_seeds2))
X <- X + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
X <- X + geom_line(data = newdata2, aes(x = avg_soil_moisture_7.6cm, 
                                       y = num_seeds2), size = 1, colour = "black")
X <- X + geom_line(data = newdata3, aes(x = avg_soil_moisture_7.6cm, 
                                        y = num_seeds2), size = 1, colour = "black", linetype = "dashed")


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
X <- X + geom_ribbon(data = newdata2, aes(ymin = LL, ymax = UL), alpha = .25)
X <- X + geom_ribbon(data = newdata3, aes(ymin = LL, ymax = UL), alpha = .25)
X


# edit data frame further for deep soil moisture. 
fitness_and_abundances_PLCO$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances_PLCO$avg_soil_moisture_12cm)
class(fitness_and_abundances_PLCO$avg_soil_moisture_12cm) 
fitness_and_abundances_VALO$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances_VALO$avg_soil_moisture_12cm)
class(fitness_and_abundances_VALO$avg_soil_moisture_12cm) 
fitness_and_abundances$avg_soil_moisture_12cm <- unfactor(fitness_and_abundances$avg_soil_moisture_12cm)
class(fitness_and_abundances$avg_soil_moisture_12cm) 

m9 <- glm.nb(num_seeds2 ~ avg_soil_moisture_12cm, 
          data = fitness_and_abundances_PLCO,  na.action=na.omit) 
summary(m9) # (avg_soil_moisture_7.6cm) is not significant
## intercept is significant = exp(3.417535)
anova(m9)
m9.1 <- lme(num_seeds2 ~ 1, 
            data = fitness_and_abundances_PLCO, random=~1|Site,  na.action=na.omit) 
summary(m9.1) # (avg_soil_moisture_7.6cm) is not significant
## intercept is significant = exp(3.61964)

m9.2 <- lme(num_seeds2 ~ avg_soil_moisture_12cm, 
            data = fitness_and_abundances_VALO2, random=~1|Site,  na.action=na.omit) 
summary(m9.2) # (soil.moisture.1.7.6cm) is not significant
## intercept is significant = exp(3.488580)
anova(m9.2)
m9.3 <- lme(num_seeds2 ~ 1, 
            data = fitness_and_abundances_VALO2, random=~1|Site,  na.action=na.omit) 
summary(m9.3) # (soil.moisture.1.7.6cm) is not significant
## intercept is significant = exp(3.54669)

newdat <- data.frame(avg_soil_moisture_12cm <- seq(15, 32, .25))
newdat$y <- exp(3.61964) 
newdat$Species = "Plectritis"

newdat2 <- data.frame(avg_soil_moisture_12cm <- seq(15, 32, .25))
newdat2$y <- exp(3.54669) 
newdat2$Species = "Valerianella"


Y <- ggplot(data = fitness_and_abundances, aes(x = as.numeric(avg_soil_moisture_12cm),
                                               y = num_seeds2))
Y <- Y + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
Y <- Y + geom_line(data = newdat, aes(x = avg_soil_moisture_12cm, y = y), size = 1)
Y <- Y + geom_line(data = newdat2, aes(x = avg_soil_moisture_12cm, y = y), size = 1, linetype = "dashed")
Y <- Y + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="Soil moisture at 12 cm depth (%VWC)", y="Predicted seeds per plant") +
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
Y