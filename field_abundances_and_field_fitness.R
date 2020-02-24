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


#################
## try P again (fitness against plectritis) with mixed effects glm.nb (glmer.nb)
str(fitness_and_abundances_PLCO2_filtered)

mixed_m1 <- glmer.nb(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants. 
                     + (1|Transect/Plot), 
          data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m1) # singular boundary due to small capture ofvarianve by random effects
# Plectritis is non-significant
## intercept is significant p < 0.001
anova(mixed_m1)

mixed_m2 <- glmer.nb(num_seeds2 ~ (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m2) # # Plectritis is non-significant


mixed_m3 <- glmer.nb(num_seeds2 ~ X.Plectritis..including.pollinator.focal.plants. 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_VALO2_filtered, na.action=na.omit) 
summary(mixed_m3)
anova(mixed_m3)


################ fitness against VALO 0.1m^2 scale.
mixed_m5 <- glmer.nb(num_seeds2 ~ X.Valerianella 
                     + (1|Transect/Plot), 
                     data = fitness_and_abundances_PLCO2_filtered,  na.action=na.omit) 
summary(mixed_m5) # singular boundary due to small capture ofvarianve by random effects

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


#############################################
# analyze relationship between species abundances
df <- read.csv("plec-and-valerianella-abundances-and-env-conditions.csv")

View(df)
df[df == "-"] <- NA
df[df == ""] <- NA
df[df == "#DIV/0!"] <- NA


# first filter out Thetis lake data
df_filtered <- df %>% 
  filter(Site != "Thetis")

# at 0.1 m^2 density

mixed_m9 <- glmer.nb(X.Valerianella ~ X.Plectritis..including.pollinator.focal.plants. 
                     + (1|Site/Transect/Plot), 
                     data = df_filtered,  na.action=na.omit) 
summary(mixed_m9) 

summary(m_abundances_1 <- glm.nb(X.Valerianella
                                 ~ X.Plectritis..including.pollinator.focal.plants., 
                                 data = df_filtered))

##########
##########

# at 1m^2 density

mixed_m10 <- glmer.nb(X.Valerianella.1m.2..including.subplots. ~ X.Plectritis.1m.2..including.subplots. 
                     + (1|Site/Transect/Plot), 
                     data = df_filtered,  na.action=na.omit) 
summary(mixed_m10) 

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
  labs(x=expression(paste("Plectritis / 1", m^2)), y=expression(paste("Valerianella / 1", m^2)))
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
S <- S + coord_fixed(ratio = 2.5, xlim = NULL, ylim = NULL, expand = TRUE)
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

ggarrange(P + theme(plot.margin = margin(r = 0, l = 0, t = 40, b = 0)), 
          T + theme(plot.margin = margin(r = 0, l = 0, t = 40, b = 0)), 
          R + theme(plot.margin = margin(r = 0, l = 0, t = 20, b = 0)), 
          S + theme(plot.margin = margin(r = 0, l = 0, t = 30, b = 0)), 
          labels = c("(A)", "(B)", "(C)", "(D)"), font.label = list(size = 16))

require(gridExtra)
PTRS <- grid.arrange(P, T, R, S, ncol=2, nrow=2, 
                   top=textGrob("", gp=gpar(fontsize=20,font=7))
)


