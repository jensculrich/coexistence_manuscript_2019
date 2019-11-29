library(ggplot2)
library(propagate)

df <- read.csv("parameter_means_sd.csv")

# P = niche differences
P <- ggplot(df, aes(x=X, y=mean_p, colour=X)) + 
  geom_errorbar(aes(ymin=mean_p-sd_p, ymax=mean_p+sd_p), width=.2, size = 1) 
P <- P + geom_point(size = 6) 
P <- P + theme_bw() 
P <- P + theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="treatment", y="niche overlap") + 
  ggtitle("(A)") 
P <- P + theme(legend.position="none")
P <- P + scale_color_manual(values=c("black", "grey"))
P <- P + theme(axis.line = element_line(size = 2))
P <- P + theme(axis.title.x = element_blank(),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 14))
P <- P + theme(plot.title = element_text(size = 20, face = "bold"))
P <- P + scale_y_continuous(limits = c(0, 1.2), breaks = c(0, .2, .4, .6, .8, 1, 1.2))
P <- P + geom_hline(yintercept= 1, linetype="dashed", 
                color = "black", size=2)
P


# Q = demographic ratio
Q <- ggplot(df, aes(x=X, y=mean_dem, colour=X)) + 
  geom_errorbar(aes(ymin=mean_dem-sd_dem, ymax=mean_dem+sd_dem), width=.2, size = 1) +
  geom_point(size = 6) 
Q <- Q + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="treatment", y="demographic ratio") + 
  ggtitle("(B)") +
  theme(legend.position="none")
Q <- Q + scale_color_manual(values=c("black", "grey"))
Q <- Q + theme(axis.line = element_line(size = 2))
Q <- Q + theme(axis.title.x = element_blank(),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 14))
Q <- Q + theme(plot.title = element_text(size = 20, face = "bold"))
Q <- Q + scale_y_continuous(limits = c(0, 2.5), breaks = c(0, .5, 1, 1.5, 2, 2.5))
Q <- Q + geom_hline(yintercept= 1, linetype="dashed", 
                    color = "black", size=2)
Q

# R = competitive response ratio
R <- ggplot(df, aes(x=X, y=mean_comp, colour=X)) + 
  geom_errorbar(aes(ymin=mean_comp-sd_comp, ymax=mean_comp+sd_comp), width=.2, size = 1) +
  geom_point(size = 6) 
R <- R + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="treatment", y="competitive response ratio") + 
  ggtitle("(C)") +
  theme(legend.position="none")
R <- R + scale_color_manual(values=c("black", "grey"))
R <- R + theme(axis.line = element_line(size = 2))
R <- R + theme(axis.title.x = element_blank(),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 14))
R <- R + theme(plot.title = element_text(size = 20, face = "bold"))
R <- R + scale_y_continuous(limits = c(0, 1.2), breaks = c(0, .2, .4, .6, .8, 1, 1.2))
R <- R + geom_hline(yintercept= 1, linetype="dashed", 
                    color = "black", size=2)
R

S <- ggplot(df, aes(x=X, y=mean_k, colour=X)) + 
  geom_errorbar(aes(ymin=mean_k-sd_k, ymax=mean_k+sd_k), width=.2, size = 1) +
  geom_point(size = 6) 
S <- S + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="treatment", y="average fitness ratio") + 
  ggtitle("(D)") +
  theme(legend.position="none")
S <- S + scale_color_manual(values=c("black", "grey"))
S <- S + theme(axis.line = element_line(size = 2))
S <- S + theme(axis.title.x = element_blank(),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 14))
S <- S + theme(plot.title = element_text(size = 20, face = "bold"))
S <- S + scale_y_continuous(limits = c(0, 1.2), breaks = c(0, .2, .4, .6, .8, 1, 1.2))
S <- S + geom_hline(yintercept= 1, linetype="dashed", 
                    color = "black", size=2)
S

library(gridExtra)
library(grid)
require(gridExtra)
grid.arrange(P, Q, R, S, ncol=2, nrow=2, 
             top=textGrob("", gp=gpar(fontsize=20,font=7))
)



#############################
#############################

fun.1 <- function(x) 1 - x
fun.2 <- function(x) 1 + x


T <- ggplot(df, aes(x=(1-mean_p), y=mean_k, colour=X)) + 
  geom_errorbar(aes(ymin=mean_k-sd_k, ymax=mean_k+sd_k), width=.05, size = 1) +
  geom_point(size = 6) 
T <- T + geom_errorbarh(aes(xmin=0, xmax=(1-mean_p+sd_p), height=.05)) 

T <- T + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x=substitute( paste("niche difference (1 - ",  italic('p'), ")" )), y="average fitness ratio") + 
  ggtitle("(A)") +
  theme(legend.position="none")
# color black for dry, grey for wet
T <- T + scale_color_manual(values=c("black", "grey"))
T <- T + theme(axis.line = element_line(size = 2))
T <- T + theme(axis.title.x = element_text(size = 20),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 16,
                 angle = 0,
                 vjust = 1
               ),
               # Y axis text
               axis.text.y = element_text(size = 16))
T <- T + theme(plot.title = element_text(size = 20, face = "bold"))
T <- T + scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = c(0, .5, 1, 1.5, 2.0))
T <- T + scale_x_continuous(expand = c(0, 0), limits = c(0, 1.05), breaks = c(0, .2, .4, .6, .8, 1))
T <- T + geom_hline(yintercept= 1, linetype="dashed", 
                    color = "black", size=1)

T <- T + stat_function(fun = fun.1)
T <- T + stat_function(fun = fun.2)


grob <- grobTree(textGrob(substitute(paste(italic("Valerianella "), "dominant")), x=0.05,  y=0.83, hjust=0,
                          gp=gpar(col="black", fontsize=20)))
grob2 <- grobTree(textGrob(substitute(paste(italic("Plectritis "), "dominant")), x=0.05,  y=0.1, hjust=0,
                          gp=gpar(col="black", fontsize=20)))
grob3 <- grobTree(textGrob("exclusion", x=0.5,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=20)))
grob4 <- grobTree(textGrob("coexistence", x=0.75,  y=0.6, hjust=0,
                          gp=gpar(col="black", fontsize=20)))
T <- T + annotation_custom(grob)
T <- T + annotation_custom(grob2)
T <- T + annotation_custom(grob3)
T <- T + annotation_custom(grob4)
T <- T + coord_fixed(ratio = 0.5)
T

############################
############################

# invasion growth rates

# wet conditions
# species j = Valerianella  
# species i = Plectritis

lambda_p <- c(428.1, 60.7)
lambda_v <- c(303.1, 58.2)
alpha_pv <- c(0.33, 0.10)
alpha_vp <- c(0.37, 0.10)
alpha_pp <- c(0.40, 0.10)
alpha_vv <- c(0.48, 0.15)

## From Hart et al. 2019
## "we quantified the invasion growth rate of sic (j) invading sic (i).
## using the following equation": 
## j_inv <- lambdaj / (1 + aji*((lambdai - 1)/aii))
## "where j is at vanishingly small density and
## N(i),t has been replaced by the expression
## which quantifies the equilibrium population
## density of i." 

## and likewise for species i
## i_inv <- lambdai / (1 + aij*((lambdaj - 1)/ajj))
## 

## invasion growth rate Valerianella (wet)
EXPR1 <- expression(lambda_v / (1 + alpha_vp*((lambda_p - 1)/alpha_pp)))
DF1 <- cbind(lambda_v, lambda_p,
             alpha_vp, alpha_pp)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat",
                  do.sim = TRUE, verbose = TRUE,
                  nsim = 1000000)
RES1

## invasion growth rate Plectritis (wet)
EXPR2 <- expression(lambda_p / (1 + alpha_pv*((lambda_v - 1)/alpha_vv)))
DF2 <- cbind(lambda_v, lambda_p,
             alpha_pv, alpha_vv)
RES2 <- propagate(expr = EXPR2, data = DF2, type = "stat",
                  do.sim = TRUE, verbose = TRUE,
                  nsim = 1000000)
RES2

# dry conditions
# species j = Valerianella  
# species i = Plectritis

# enter parameter estimates with sd's
lambda_p_dry <- c(364.2, 42.51)
lambda_v_dry <- c(636.1, 89.7)
alpha_pv_dry <- c(0.30, 0.10)
alpha_vp_dry <- c(0.43, 0.11)
alpha_pp_dry <- c(0.23, 0.05)
alpha_vv_dry <- c(0.72, 0.15)


## invasion growth rate Valerianella (dry)
EXPR1_dry <- expression(lambda_v_dry / (1 + alpha_vp_dry*((lambda_p_dry - 1)/alpha_pp_dry)))
DF1_dry <- cbind(lambda_v_dry, lambda_p_dry,
             alpha_vp_dry, alpha_pp_dry)
RES1_dry <- propagate(expr = EXPR1_dry, data = DF1_dry, type = "stat",
                  do.sim = TRUE, verbose = TRUE,
                  nsim = 1000000)
RES1_dry

## invasion growth rate Plectritis (dry)
EXPR2_dry <- expression(lambda_p_dry / (1 + alpha_pv_dry*((lambda_v_dry - 1)/alpha_vv_dry)))
DF2_dry <- cbind(lambda_v_dry, lambda_p_dry,
             alpha_pv_dry, alpha_vv_dry)
RES2_dry <- propagate(expr = EXPR2_dry, data = DF2_dry, type = "stat",
                  do.sim = TRUE, verbose = TRUE,
                  nsim = 1000000)
RES2_dry


mean_inv_growth_rate_j <- as.data.frame(rbind(0.765, 0.935))
sd_inv_growth_rate_j <- as.data.frame(rbind(0.335, 0.357))
mean_inv_growth_rate_i <- as.data.frame(rbind(2.05, 1.37))
sd_inv_growth_rate_i <- as.data.frame(rbind(1.01, 0.593))

df2 <- merge(df, mean_inv_growth_rate_j)

#############################
# plot invasion growth rates

E <- ggplot(df, aes(x=X, y=mean_p, colour=X)) + 
  geom_errorbar(aes(ymin=mean_p-sd_p, ymax=mean_p+sd_p), width=.2, size = 1) 
E <- E + geom_point(size = 6) 
E <- E + theme_bw() 
E <- E + theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="treatment", y="niche overlap") + 
  ggtitle("(A)") 
E <- E + theme(legend.position="none")
E <- E + scale_color_manual(values=c("black", "grey"))
E <- E + theme(axis.line = element_line(size = 2))
E <- E + theme(axis.title.x = element_blank(),
               # Y axis title
               axis.title.y = element_text(size = 20),
               # X axis text
               axis.text.x = element_text(
                 size = 20,
                 angle = 0,
                 vjust = .5
               ),
               # Y axis text
               axis.text.y = element_text(size = 14))
E <- E + theme(plot.title = element_text(size = 20, face = "bold"))
E <- E + scale_y_continuous(limits = c(0, 1.2), breaks = c(0, .2, .4, .6, .8, 1, 1.2))
E <- E + geom_hline(yintercept= 1, linetype="dashed", 
                    color = "black", size=2)
E