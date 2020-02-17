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


fitness_and_abundances <- read.csv("fitness_and_abundances.csv")
fitness_and_abundances_PLCO <- read.csv("fitness_and_abundances_PLCO.csv")
fitness_and_abundances_VALO <- read.csv("fitness_and_abundances_VALO.csv")

### remove the big outlier from the data set (many more seeds preditced vs any other plot)
which.max(fitness_and_abundances_PLCO$num_seeds2)
fitness_and_abundances_PLCO2 <- fitness_and_abundances_PLCO[-31, ]

### remove the big outlier from the data set (more than twice as many Plectritis density than second highest value)
which.max(fitness_and_abundances_VALO$X.Plectritis..including.pollinator.focal.plants.)
fitness_and_abundances_VALO2 <- fitness_and_abundances_VALO[-9, ]

### at the 0.1m^2 scale at least the number of seeds per plant is highly variable at low abundances
### but then regulates as the number of plectritis passes ~15 per 0.1m^2
### maybe suggesting that population regulates itself as densities increaseS
fitness_and_abundances_VALO2_filtered <- fitness_and_abundances_VALO2 %>% 
  filter(Site != "Thetis")
fitness_and_abundances_PLCO2_filtered <- fitness_and_abundances_PLCO2 %>% 
  filter(Site != "Thetis")

m1 <- lme(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants., 
          data = fitness_and_abundances_PLCO2_filtered, random=~1|Site,  na.action=na.omit) 
summary(m1) # # Plectritis is non-significant
## intercept is significant is exp(3.5716245) = 35.56037, p < 0.001
m1.2 <- lme(num_seeds2 ~ 1, 
            data = fitness_and_abundances_PLCO2_filtered, random=~1|Site,  na.action=na.omit) 
# reduce the model to just calc intercept.
summary(m1.2) # 
## intercept is significant = 37.3240, p < 0.001, SSE = 6.77.
anova(m1)

m1.3 <- lme(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants., 
            data = fitness_and_abundances_VALO2_filtered, random=~1|Site,  na.action=na.omit) 
m1.4 <- lme(num_seeds2 ~ 1, 
            data = fitness_and_abundances_VALO2_filtered, random=~1|Site,  na.action=na.omit) 
summary(m1.3) # Plectritis is significant = -0.20090, p <  0.00299 
# intercept is significant = 36.80011, p < 0.05
# R^2 is low, 0.1059, SSE = 6.626
anova(m1.3)
resid(m1.3)

# calculate residuals for PLCO
fitness_and_abundances_PLCO2_filtered <- fitness_and_abundances_PLCO2 %>%
  mutate(residual = num_seeds2 - 37.3240) 
# calculate residuals for VALO
r <- fitness_and_abundances_VALO2_filtered %>%
  do(augment(lme(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants., 
                 data = ., random=~1|Site,  na.action=na.omit)))

fitness_and_abundances_VALO2_filtered <- fitness_and_abundances_VALO2_filtered %>%
  mutate(residual = resid(m1.3)) 


# calculate some new data to show relationship between VALO fitness and PLCO abundance
newdat <- data.frame(X.Plectritis..including.pollinator.focal.plants. <- seq(0, 60, .5))
newdat <- newdat %>% 
  mutate(y = 36.80011 - 0.21019*X.Plectritis..including.pollinator.focal.plants.) %>% 
  mutate(Species = "VALO")

fitness_and_abundances_filtered <- fitness_and_abundances %>% 
  filter(Site != "Thetis")

levels(fitness_and_abundances_filtered$Species)[levels(fitness_and_abundances_filtered$Species)=="PLCO"] <- "Plectritis"
levels(fitness_and_abundances_filtered$Species)[levels(fitness_and_abundances_filtered$Species)=="VALO"] <- "Valerianella"

P <- ggplot(data = fitness_and_abundances_filtered, aes(x = X.Plectritis..including.pollinator.focal.plants.,
                                                        y = num_seeds2), color = Species)
P <- P + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
P <- P + geom_segment(aes(x=0,xend=130,y=37.3240,yend=37.3240), color ="black", size = 1)
P <- P + geom_line(data = newdat, aes(x = X.Plectritis..including.pollinator.focal.plants., y = y, color = Species), size = 1, linetype = "dashed")
P <- P + theme_bw() + 
  theme(legend.position = c(0.75, 0.9), plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  labs(x=expression(paste("Plectritis density (0.1", m^2, ")")), y="Predicted seeds per plant") 
P <- P + scale_color_manual(values=c("black", "black"))
P <- P + scale_shape_manual(values=c(19, 1))
P <- P + theme(axis.line = element_line(size = 2))
P <- P + theme(legend.title=element_blank())
P <- P + theme(legend.text=element_text(size=rel(2)))
# X axis title
P <- P + theme(axis.title.x = element_text(vjust = 0,
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
P <- P  + theme(plot.title = element_text(size = 20, face = "bold"))
P <- P + guides(color = FALSE)
P <- P + theme(axis.line = element_line(colour = 'black', size = 1))
P <- P + theme(plot.margin = margin(r = 4, l = 4))
P


################ fitness against VALO 0.1m^2 scale.
m5 <- lme(num_seeds2 ~  X.Valerianella, 
          data = fitness_and_abundances_PLCO2_filtered, random=~1|Site,  na.action=na.omit) 
summary(m5) # # Valerianella is non-significant
## intercept is significant = 36.06412, p < 0.001

m5.2 <- lme(num_seeds2 ~ 1, 
            data = fitness_and_abundances_PLCO2_filtered, random=~1|Site,  na.action=na.omit) 
summary(m5.2) # # Valerianella removed
## intercept is significant = 37.32403, p < 0.001, SSE = 6.77
anova(m5)

m5.3 <- lme(num_seeds2 ~ X.Valerianella, 
            data = fitness_and_abundances_VALO2_filtered, random=~1|Site,  na.action=na.omit) 
summary(m5.3) # Valerianella is significant = 0.21199, p < 0.05
# intercept is significant = 35.45480, p < 0.001
# R^2 is low, 0.07606, SSE = 6.681
anova(m5.3)

# calculate some new data to sshow relationship between VALO fitness and VALO abundance
newdat2 <- data.frame(X.Valerianella <- seq(0, 46.5, .5))
newdat2 <- newdat2 %>% 
  mutate(y = 35.45480 + 0.21199*X.Valerianella) %>% 
  mutate(Species = "VALO")

# calculate residuals
fitness_and_abundances_PLCO2_filtered <- fitness_and_abundances_PLCO2_filtered %>%
  mutate(residual = num_seeds2 - 37.32403) 

T <- ggplot(data = fitness_and_abundances_filtered, aes(x = X.Valerianella,
                                                        y = num_seeds2), color = Species)
T <- T + geom_point(size = 2, alpha = 0.7, aes(shape = Species))
T <- T + geom_segment(aes(x=0,xend=20,y=37.32403,yend=37.32403), color ="black", size = 1)
T <- T + geom_line(data = newdat2, aes(x = X.Valerianella, y = y, color = Species), size = 1, linetype = "dashed")
T <- T + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x=expression(paste("Valerianella density (0.1", m^2, ")")), y="Predicted seeds per plant") +
  theme(legend.position="none")
T <- T + scale_shape_manual(values=c(19, 1))
T <- T + scale_color_manual(values=c("black", "black"))
T <- T + theme(axis.line = element_line(size = 2))
T <- T + theme(axis.title.x = element_text(vjust = 0,
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
T <- T + theme(plot.title = element_text(size = 20, face = "bold"))
T <- T + theme(axis.line = element_line(colour = 'black', size = 1))
T

m6 <- lme(log(abs(residual)+1) ~ log(X.Valerianella+1), 
          data = fitness_and_abundances_PLCO2_filtered, random=~1|Site,  na.action=na.omit) 
summary(m6) # only intercept is significant
# r^2 is VERY low 
anova(m6)

require(gridExtra)
PT <- grid.arrange(P, T, ncol=1, nrow=2, 
                   top=textGrob("", gp=gpar(fontsize=20,font=7))
)


# analyze relationship between species abundances
df <- read.csv("plec-and-valerianella-abundances-and-env-conditions.csv")

View(df)
df[df == "-"] <- NA
df[df == ""] <- NA
df[df == "#DIV/0!"] <- NA


# first filter out Thetis lake data
df_filtered <- df %>% 
  filter(Site != "Thetis")
# fit a non-linear  function to describe relationship 
#nLL <- function(beta1, beta2, beta3) -sum(stats::dpois(
#            df_filtered$X.Valerianella/(beta1+beta2*df_filtered$X.Plectritis..including.pollinator.focal.plants.)+beta3*df_filtered$X.Plectritis..including.pollinator.focal.plants.^2), log = TRUE)
#mle_fit <- mle(nLL, start = list(beta1 = 1, beta2 = 1, beta3 = 1), nobs = NROW(df_filtered))


newdat3 <- data.frame(X.Plectritis..including.pollinator.focal.plants. <- seq(0, 4.5, 0.1))
newdat3$y <- 1.78444 - 0.37723*newdat3$X.Plectritis..including.pollinator.focal.plants.....seq.0..4.5.. 

R <- ggplot(df_filtered, aes(x = log(1 + X.Plectritis..including.pollinator.focal.plants.), 
                             y = log(1 + X.Valerianella)))
R <- R + geom_point(aes(log(1 + X.Plectritis..including.pollinator.focal.plants.), 
                        log(1 + X.Valerianella), shape = Site), size = 2, alpha = 1) 
R <- R + geom_line(data = newdat3, aes(x = X.Plectritis..including.pollinator.focal.plants.....seq.0..4.5.., 
                                       y = y), size = 1, colour = "black")
R <- R + theme_bw() + theme(plot.title = element_text(hjust = 0), 
                            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.position = c(.8, .99), legend.title=element_text()) + 
  labs(x=expression(paste("log(Plectritis density per 0.1", m^2, ")")), y=expression(paste("log(Valerianella density per 0.1", m^2, ")")))
R <- R + scale_color_manual(values=c("black", "black"))
R <- R + scale_shape_manual(values=c(1, 16))
R <- R + theme(axis.line = element_line(size = 2))
R <- R + theme(axis.title.x = element_text(vjust = 0,
                                           size = 18),
               # Y axis title
               axis.title.y = element_text(size = 18),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 20))
R <- R + theme(plot.title = element_text(size = 20, face = "bold"))
R <- R + theme(legend.text=element_text(size=rel(2)))
R <- R + theme(legend.title=element_text(size=rel(2)))
R <- R + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
R <- R + theme(axis.line = element_line(colour = 'black', size = 1))
R

m2 <- lm(log(1 + X.Valerianella.1m.2..including.subplots.) ~ log(1 + X.Plectritis.1m.2..including.subplots.)
         , data = df_filtered)
summary(m2)
m2.1 <- lm(log(1 + X.Valerianella.1m.2..including.subplots.) ~ log(1 + X.Plectritis.1m.2..including.subplots.) + as.numeric(avg_soil_moisture_7.6cm)
           , data = df_filtered)
summary(m2.1)

newdat4 <- data.frame(X.Plectritis.1m.2..including.subplots. <- seq(0, 7, 0.1))
newdat4$y <- 3.7998 - 0.2947*newdat4$X.Plectritis.1m.2..including.subplots.....seq.0..7..0.1. 

S <- ggplot(df_filtered, aes(x = log(1 + X.Plectritis.1m.2..including.subplots.), 
                             y = log(1 + X.Valerianella.1m.2..including.subplots.)))
S <- S + geom_point(aes(log(1 + X.Plectritis.1m.2..including.subplots.), 
                        log(1 + X.Valerianella.1m.2..including.subplots.), shape = Site), size = 2, alpha = 0.7) 
S <- S + geom_line(data = newdat4, aes(x = X.Plectritis.1m.2..including.subplots.....seq.0..7..0.1., 
                                       y = y), size = 1, colour = "black")
S <- S + theme_bw() + theme(plot.title = element_text(hjust = 0), 
                            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.position = "none", legend.title=element_blank()) + 
  labs(x=expression(paste("log(Plectritis density per 1", m^2, ")")), y=expression(paste("log(Valerianella density per 1", m^2, ")")))
S <- S + scale_color_manual(values=c("black", "black"))
S <- S + scale_shape_manual(values=c(1, 16))
S <- S + theme(axis.line = element_line(size = 2))
S <- S + theme(axis.title.x = element_text(vjust = 0,
                                           size = 18),
               # Y axis title
               axis.title.y = element_text(size = 18),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 20))
S <- S + theme(plot.title = element_text(size = 20, face = "bold"))
S <- S + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
S <- S + theme(axis.line = element_line(colour = 'black', size = 1))
S


# reload figure P first
ggarrange(P + theme(plot.margin = margin(r = 0, l = 10, t = 80, b = 20)), 
          T + theme(plot.margin = margin(r = 0, l = 10, t = 80, b = 20)), 
          R + theme(plot.margin = margin(r = 0, l = 10, t = 20, b = 0)), 
          S + theme(plot.margin = margin(r = 0, l = 30, t = 20, b = 0)), 
          labels = c("(A)", "(B)", "(C)", "(D)"), font.label = list(size = 20))


