# this file analyzes relationships between fecundity 
# and hetero-/con-specific densities
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


fitness_and_abundances_filtered <- fitness_and_abundances %>% 
  filter(Site != "Thetis")

levels(fitness_and_abundances_filtered$Species)[levels(fitness_and_abundances_filtered$Species)=="PLCO"] <- "Plectritis"
levels(fitness_and_abundances_filtered$Species)[levels(fitness_and_abundances_filtered$Species)=="VALO"] <- "Valerianella"



#################
## try P again (fitness against plectritis) with mixed effects glm.nb (glmer.nb)
str(fitness_and_abundances_PLCO2_filtered)
mixed_m1 <- glmer.nb(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants. 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m1) # # Plectritis is non-significant
## intercept is significant is 36.22965, p < 0.001

summary(mixed_m1) # singular boundary due to small capture ofvarianve by random effects
# Plectritis is non-significant
## intercept is significant p < 0.001
anova(mixed_m1)

mixed_m2 <- glmer.nb(num_seeds2 ~ (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m2) # # Plectritis is non-significant
# reduce the model to just calc intercept.
summary(m1.2) # 
## intercept is significant = 37.35814, p < 0.001, SSE = 6.77.
anova(m1)


mixed_m3 <- glm.nb(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants., 
                     data = fitness_and_abundances_VALO2_filtered,  na.action=na.omit) 
summary(mixed_m3)
anova(mixed_m3)


#only need to graph negative valerianella because no relationship with plectritis fitness
newdata <- data.frame(
  X.Plectritis..including.pollinator.focal.plants. = 
    rep(seq(from = min(fitness_and_abundances_VALO2_filtered$X.Plectritis..including.pollinator.focal.plants.), 
            to = 60, length.out = 100)))
newdata <- cbind(newdata, predict(mixed_m3, newdata, type = "link", se.fit=TRUE))
newdata <- within(newdata, {
  num_seeds2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})



P <- ggplot(fitness_and_abundances_filtered, aes(x = X.Plectritis..including.pollinator.focal.plants., 
                             y = num_seeds2))
P <- P + geom_point(aes(X.Plectritis..including.pollinator.focal.plants., 
                        num_seeds2, shape = Species), size = 2, alpha = 0.7) 
P <- P + geom_line(data = newdata, aes(x = X.Plectritis..including.pollinator.focal.plants., 
                                        y = num_seeds2), size = 1, colour = "black", linetype = "dashed")
P <- P + geom_segment(aes(x=0,xend=110,y=37.1,yend=37.1), color ="black", size = 1)
P <- P + theme_bw() + theme(plot.title = element_text(hjust = 0), 
                            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.position = "none", legend.title=element_text()) + 
  labs(x=expression(paste("P. congesta / 0.1", m^2)), y="Seeds / plant")
P <- P + scale_color_manual(values=c("black", "black"))
P <- P + scale_shape_manual(values=c(19, 1))
P <- P + theme(axis.line = element_line(size = 2))
P <- P + theme(axis.title.x = element_text(vjust = 0,
                                           size = 16),
               # Y axis title
               axis.title.y = element_text(size = 16),
               # X axis text
               axis.text.x = element_text(
                 size = 16,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 16))
P <- P + theme(plot.title = element_text(size = 20, face = "bold"))
P <- P + theme(axis.line = element_line(colour = 'black', size = 1))
P <- P + geom_ribbon(data = newdata, aes(ymin = LL, ymax = UL), alpha = .25)
P <- P + theme(legend.text=element_text(size=rel(1.5)))
P <- P + theme(legend.title=element_text(size=rel(1.5)))
P <- P + ggtitle("(A)")
P


################ fitness against VALO 0.1m^2 scale.

mixed_m5 <- glmer.nb(num_seeds2 ~ X.Valerianella 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m5) # singular boundary due to small capture ofvarianve by random effects

m6 <- lme(log(abs(residual)+1) ~ log(X.Valerianella+1), 
          data = fitness_and_abundances_PLCO2_filtered, random=~1|Site,  na.action=na.omit) 
summary(m6) # only intercept is significant
# r^2 is VERY low 

anova(m6)
mixed_m6 <- glmer.nb(num_seeds2 ~ (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m6) # singular boundary due to small capture ofvarianve by random effects


mixed_m7 <- glmer.nb(num_seeds2 ~ X.Valerianella 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_VALO2_filtered,  na.action=na.omit) 
summary(mixed_m7) 

mixed_m8 <- glmer.nb(num_seeds2 ~ (1|Transect/Plot), 
                     data = fitness_and_abundances_VALO2_filtered,  na.action=na.omit) 
summary(mixed_m8) 

T <- ggplot(data = fitness_and_abundances_filtered, aes(x = X.Valerianella,
                                                        y = num_seeds2), color = Species)
T <- T + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
T <- T + geom_segment(aes(x=0,xend=50,y=37.1,yend=37.1), color ="black", size = 1)
T <- T + geom_segment(aes(x=0,xend=50,y=34.3,yend=34.3), color = "black", linetype = "dashed", size = 1)
T <- T + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x=expression(paste("V. locusta / (0.1", m^2, ")")), y="Seeds / plant") +
  theme(legend.position="none")
T <- T + scale_shape_manual(values=c(19, 1))
T <- T + scale_color_manual(values=c("black", "black"))
T <- T + theme(axis.line = element_line(size = 2))
T <- T + theme(axis.title.x = element_text(vjust = 0,
                                           size = 16),
               # Y axis title
               axis.title.y = element_text(size = 16),
               # X axis text
               axis.text.x = element_text(
                 size = 16,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 16))
T <- T + theme(plot.title = element_text(size = 20, face = "bold"))
T <- T + theme(axis.line = element_line(colour = 'black', size = 1))
T <- T + ggtitle("(B)")
T



#############################################
# analyze relationship between species abundances
df <- read.csv("plec-and-valerianella-abundances-and-env-conditions.csv")

df[df == "-"] <- NA
df[df == ""] <- NA
df[df == "#DIV/0!"] <- NA

df_filtered <- df %>%
  filter(Site != "Thetis")
  # at 0.1 m^2 density
  
mixed_m9 <- glmer.nb(X.Valerianella ~ X.Plectritis..including.pollinator.focal.plants. 
                       + (1|Site/Transect/Plot), 
                       data = df_filtered,  na.action=na.omit) 
summary(mixed_m9) 
anova(mixed_m9)

summary(m_abundances_1 <- glm.nb(X.Valerianella
                                 ~ X.Plectritis..including.pollinator.focal.plants., 
                                 data = df_filtered))

newdata3 <- data.frame(
  X.Plectritis..including.pollinator.focal.plants. = 
    rep(seq(from = min(df_filtered$X.Plectritis..including.pollinator.focal.plants.), 
            to = max(df_filtered$X.Plectritis..including.pollinator.focal.plants.), length.out = 100)))

newdata3 <- cbind(newdata3, predict(m_abundances_1, newdata3, type = "link", se.fit=TRUE))

newdata3 <- within(newdata3, {
  X.Valerianella <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


R <- ggplot(df_filtered, aes(x = X.Plectritis..including.pollinator.focal.plants., 
                             y = X.Valerianella))
R <- R + geom_point(data = df_filtered, aes(X.Plectritis..including.pollinator.focal.plants., 
                                            X.Valerianella, shape = Site), size = 2, alpha = 1) 
R <- R + geom_line(data = newdata3, aes(x = X.Plectritis..including.pollinator.focal.plants., 
                                        y = X.Valerianella), size = 1, colour = "black")
R <- R + theme_bw() + theme(plot.title = element_text(hjust = 0), 
                            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.position = c(.66, .7), legend.title=element_text()) + 
  labs(x=expression(paste("P. congesta / 0.1", m^2)), y=expression(paste("V. locusta / 0.1", m^2)))
R <- R + scale_color_manual(values=c("black", "black"))
R <- R + scale_shape_manual(values=c(1, 16))
R <- R + theme(axis.line = element_line(size = 2))
R <- R + theme(axis.title.x = element_text(vjust = 0,
                                           size = 16),
               # Y axis title
               axis.title.y = element_text(size = 16),
               # X axis text
               axis.text.x = element_text(
                 size = 16,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 16))
R <- R + theme(plot.title = element_text(size = 20, face = "bold"))
R <- R + theme(legend.text=element_text(size=rel(1.5)))
R <- R + theme(legend.title=element_text(size=rel(1.5)))
R <- R + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
R <- R + theme(axis.line = element_line(colour = 'black', size = 1))
R <- R + geom_ribbon(data = newdata3, aes(ymin = LL, ymax = UL), alpha = .25)
R

##########
##########

# at 1m^2 density

mixed_m10 <- glmer.nb(X.Valerianella.1m.2..including.subplots. ~ X.Plectritis.1m.2..including.subplots. 
                      + (1|Site/Transect/Plot), 
                      data = df_filtered,  na.action=na.omit) 
summary(mixed_m10) 
anova(mixed_m10)

summary(m_abundances_2 <- glm.nb(X.Valerianella.1m.2..including.subplots.
                                 ~ X.Plectritis.1m.2..including.subplots., 
                                 data = df_filtered))
newdata4 <- data.frame(
  X.Plectritis.1m.2..including.subplots. = 
    rep(seq(from = min(df_filtered$X.Plectritis.1m.2..including.subplots.), 
            to = max(df_filtered$X.Plectritis.1m.2..including.subplots.), length.out = 100)))
newdata4 <- cbind(newdata4, predict(m_abundances_2, newdata4, type = "link", se.fit=TRUE))
newdata4 <- within(newdata4, {
  X.Valerianella.1m.2..including.subplots. <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata4, aes(X.Plectritis.1m.2..including.subplots., 
                     X.Valerianella.1m.2..including.subplots.)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line(aes(), size = 2) +
  labs(x = "x", y = "y")

S <- ggplot(df_filtered, aes(x = X.Plectritis.1m.2..including.subplots., 
                             y = X.Valerianella.1m.2..including.subplots.))
S <- S + geom_point(aes(X.Plectritis.1m.2..including.subplots., 
                        X.Valerianella.1m.2..including.subplots., shape = Site), size = 2, alpha = 0.7) 
S <- S + geom_line(data = newdata4, aes(x = X.Plectritis.1m.2..including.subplots., 
                                        y = X.Valerianella.1m.2..including.subplots.), size = 1, colour = "black")
S <- S + theme_bw() + theme(plot.title = element_text(hjust = 0), 
                            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.position = "none", legend.title=element_blank()) + 
  labs(x=expression(paste("P. congesta / 1", m^2)), y=expression(paste("V. locusta / 1", m^2)))
S <- S + scale_color_manual(values=c("black", "black"))
S <- S + scale_shape_manual(values=c(1, 16))
S <- S + theme(axis.line = element_line(size = 2))
S <- S + theme(axis.title.x = element_text(vjust = 0,
                                           size = 16),
               # Y axis title
               axis.title.y = element_text(size = 16),
               # X axis text
               axis.text.x = element_text(
                 size = 16,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 16))
S <- S + theme(plot.title = element_text(size = 20, face = "bold"))
S <- S + coord_fixed(ratio = 1.5, xlim = NULL, ylim = NULL, expand = TRUE)
S <- S + theme(axis.line = element_line(colour = 'black', size = 1))
S <- S + geom_ribbon(data = newdata4, aes(ymin = LL, ymax = UL), alpha = .25)
S

############
# reload figure P first
ggarrange(P + theme(plot.margin = margin(r = 0, l = 35, t = 25, b = 0)), 
          T + theme(plot.margin = margin(r = 2, l = 35, t = 25, b = 0)), 
          R + theme(plot.margin = margin(r = 0, l = 10, t = 20, b = 10)), 
          S + theme(plot.margin = margin(r = 2, l = 10, t = 20, b = 10)), 
          labels = c("(A)", "(B)", "(C)", "(D)"), font.label = list(size = 16))
ggarrange(P + theme(plot.margin = margin(r = 5, l = 5, t = 5, b = 5)), 
          T + theme(plot.margin = margin(r = 5, l = 5, t = 5, b = 5)), 
          R + theme(plot.margin = margin(r = 5, l = 5, t = 5, b = 5)), 
          S + theme(plot.margin = margin(r = 5, l = 5, t = 5, b = 5)), 
          labels = c("(A)", "(B)", "(C)", "(D)"), font.label = list(size = 16))
require(gridExtra)
PTRS <- grid.arrange(P, T, R, S, ncol=2, nrow=2, 
                     top=textGrob("", gp=gpar(fontsize=20,font=7))
)

grid.arrange(R, S, ncol=1, nrow=2, 
             top=textGrob("", gp=gpar(fontsize=20,font=7))
)


#####################################

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


m7 <- lme(num_seeds2 ~  as.numeric(X.grasscover.subplot), 
          data = fitness_and_abundances_PLCO,random=~1|Site,  na.action=na.omit) 
summary(m7) # X.grasscover.subplot is minorly significant (positive)
## intercept is significant =
anova(m7)
m7.1 <- lm(num_seeds2 ~  as.numeric(X.grasscover.subplot), 
           data = fitness_and_abundances_VALO) 
summary(m7.1) # X.grasscover.subplot is minorly significant (positive)
## intercept is significant =


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
X <- X + ggtitle("(C)")
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
Y <- Y + ggtitle("(D)")
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


fitness_and_abundances_temp3 <- fitness_and_abundances_temp2[!is.na(fitness_and_abundances_temp2$soil.depth.1),]
fitness_and_abundances_temp3 <- na.omit(fitness_and_abundances_temp2)
!is.na(fitness_and_abundances_temp2$soil.depth.1)

Z <- ggerrorplot(fitness_and_abundances_temp3, x = "soil.depth.1", y = "num_seeds2", 
                  desc_stat = "mean_se", color  = "Species", position = position_dodge(0.3) 
)
Z <- Z + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="soil depth", y="Seeds / plant") +
  theme(legend.position="none")
Z <- Z + scale_color_manual(values=c("black", "grey"))
Z <- Z + scale_shape_manual(values=c(1, 19))
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
Z <- Z + ggtitle("(E)")

Z  
  
 
############
############


require(gridExtra)
grid.arrange(P, T, X, Y, Z, ncol=2, nrow=3, 
             top=textGrob("", gp=gpar(fontsize=20,font=7))
)

ggarrange(#V + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
  Z + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
  X + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
  Y + theme(plot.margin = margin(r = 10, l = 10, t = 60, b = 10)), 
  labels = c("(A)", "(B)", "(C)", "(D)"), font.label = list(size = 20))

