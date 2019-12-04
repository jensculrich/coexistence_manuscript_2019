# program code to analyze relationship between species abundances

library(tidyverse)
library(purrr)
library(repurrrsive)
library(listviewer)
library(stats4)
library(nlme)

df <- read.csv("plec-and-valerianella-abundances-and-env-conditions.csv")

View(df)
df[df == "-"] <- NA
df[df == ""] <- NA
df[df == "#DIV/0!"] <- NA

# when PLCO is absent or low Valerianella is very high, 
# when PLCO is high Valerianella abundances are low.
#par(mfrow=c(1,2),mai=c(0.6,0.6,0.6,0.6),cex=0.55)
#plot(X.Valerianella ~ X.Plectritis..including.pollinator.focal.plants., data = df, 
#     xlab = "# of Plectritis", ylab = "# Valerianella", main = "Abundances at 0.1m^2")
#plot(X.Valerianella.1m.2..including.subplots. ~ X.Plectritis.1m.2..including.subplots. , data = df, 
#     xlab = "# of Plectritis", ylab = "# Valerianella", main = "Abundances at 1m^2")
#lm1 <- glm(X.Valerianella ~ X.Plectritis..including.pollinator.focal.plants.
#          , data = df, family = poisson)
#summary(lm1)
#glm1 <- glm(df$X.Valerianella~df$X.Plectritis..including.pollinator.focal.plants., family=poisson()) 
#summary(glm1)

#df[df == "NA"] <- 0
#nLL <- function(abundance, alpha) -sum(stats::dpois(df$X.Valerianella, 
#      abundance/(1+alpha*df$X.Plectritis..including.pollinator.focal.plants.), log = TRUE))
# mle_fit <- mle(nLL, start = list(abundance = 30, alpha = 0.01), nobs = NROW(df))

#newdat = data.frame(Plectritis = df$X.Plectritis..including.pollinator.focal.plants.)
#predicted_values <- predict(glm1, newdata = newdat, interval="confidence")

# first filter out Thetis lake data
df_filtered <- df %>% 
  filter(Site != "Thetis")
# fit a non-linear second degree inverse polynomial function to describe relationship 
#nLL <- function(beta1, beta2, beta3) -sum(stats::dpois(
#            df_filtered$X.Valerianella/(beta1+beta2*df_filtered$X.Plectritis..including.pollinator.focal.plants.)+beta3*df_filtered$X.Plectritis..including.pollinator.focal.plants.^2), log = TRUE)
#mle_fit <- mle(nLL, start = list(beta1 = 1, beta2 = 1, beta3 = 1), nobs = NROW(df_filtered))


newdat2 <- data.frame(X.Plectritis..including.pollinator.focal.plants. <- seq(0, 4.5, 0.1))
newdat2$y <- 1.78444 - 0.37723*newdat2$X.Plectritis..including.pollinator.focal.plants.....seq.0..4.5.. 

R <- ggplot(df_filtered, aes(x = log(1 + X.Plectritis..including.pollinator.focal.plants.), 
                             y = log(1 + X.Valerianella)))
R <- R + geom_point(aes(log(1 + X.Plectritis..including.pollinator.focal.plants.), 
                        log(1 + X.Valerianella), shape = Site), size = 3, alpha = 1) 
R <- R + geom_line(data = newdat2, aes(x = X.Plectritis..including.pollinator.focal.plants.....seq.0..4.5.., 
                                       y = y), size = 1.5, colour = "black")
R <- R + theme_bw() + theme(plot.title = element_text(hjust = 0), 
                            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.position = c(.8, .99), legend.title=element_text()) + 
  labs(x=expression(paste("log(Plectritis density per 0.1", m^2, ")")), y=expression(paste("log(Valerianella density per 0.1", m^2, ")"))) + 
  ggtitle("(A)")
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
                        log(1 + X.Valerianella.1m.2..including.subplots.), shape = Site), size = 3, alpha = 0.7) 
S <- S + geom_line(data = newdat4, aes(x = X.Plectritis.1m.2..including.subplots.....seq.0..7..0.1., 
                                       y = y), size = 1.5, colour = "black")
S <- S + theme_bw() + theme(plot.title = element_text(hjust = 0), 
                            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.position = "none", legend.title=element_blank()) + 
  labs(x=expression(paste("log(Plectritis density per 1", m^2, ")")), y=expression(paste("log(Valerianella density per 1", m^2, ")"))) + 
  ggtitle("(B)")
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
S

require(grid)
require(gridExtra)
grid.arrange(R, S, ncol=1, nrow=2, 
             top=textGrob("", gp=gpar(fontsize=20,font=7)))
