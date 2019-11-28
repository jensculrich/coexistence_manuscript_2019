library(ggplot2)

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

