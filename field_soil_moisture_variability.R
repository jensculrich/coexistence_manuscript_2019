########################
########################
# soil moisture spatial variation

########################
library(tidyverse)

# load dataframe for 12cm depth soil moisture readings
# and for 7cm depth soil moisture readings

data_long_12cm <- read.csv("soil_moisture_12cm.csv")
data_long_7cm <- read.csv("soil_moisture_7cm.csv")

########################
library(nlme) 
plot.lme <- lme(moisture~1, random=~1|plot_id, data=data_long_7cm, na.action=na.omit) 
VarCorr(plot.lme) 
# intercept gives between plot variation, residual gives within plot variation.

plot.lme <- lme(moisture~1, random=~1|plot_id, data=data_long_12cm, na.action=na.omit) 
VarCorr(plot.lme) 

# intercept gives between plot variation, residual gives within plot variation.
Transect.lme <- lme(moisture~1, random=~1|Transect, data=data_long_7cm, na.action=na.omit) 
VarCorr(Transect.lme) 

# nested
nested.lme_7cm <- lme(moisture~1, random=~1| Transect / plot_id, data=data_long_7cm, na.action=na.omit) 
x <- VarCorr(nested.lme_7cm)

nested.lme_12cm <- lme(moisture~1, random=~1| Transect / plot_id, data=data_long_12cm, na.action=na.omit) 
y <- VarCorr(nested.lme_12cm) 

###########
# restructure VarCorr.lme's into plotable data frames

df_7cm <- x
df_7cm <- df_7cm[-c(1, 3), ] 
rownames(df_7cm) <- c("between transects", "between plots", "within plots")
df_7cm <- as.data.frame(df_7cm)
df_7cm$Variance <- as.numeric(as.character(df_7cm$Variance))
df_7cm$StdDev <- as.numeric(as.character(df_7cm$StdDev))
df_7cm$Depth <- as.factor(c(7, 7, 7))
df_7cm$proportion_of_variance <- (df_7cm$Variance/(sum(df_7cm$Variance)))
df_7cm$proportional_stdev <- (df_7cm$StdDev/(sum(df_7cm$Variance)))
str(df_7cm)

df_12cm <- y
df_12cm <- df_12cm[-c(1, 3), ] 
rownames(df_12cm) <- c("between transects", "between plots", "within plots")
df_12cm <- as.data.frame(df_12cm)
df_12cm$Variance <- as.numeric(as.character(df_12cm$Variance))
df_12cm$StdDev <- as.numeric(as.character(df_12cm$StdDev))
df_12cm$Depth <- as.factor(c(12, 12, 12))
df_12cm$proportion_of_variance <- (df_12cm$Variance/(sum(df_12cm$Variance)))
df_12cm$proportional_stdev <- (df_12cm$StdDev/(sum(df_12cm$Variance)))
str(df_12cm)

# combine 7cm and 12cm data frames
df <- rbind(df_7cm, df_12cm)
df$Scale <- c("between transects", "between plots", "within plots", 
              "between transects", "between plots", "within plots")
rownames(df) <- c()

df


######

P <- ggplot(df, aes(x=Scale, y=proportion_of_variance, colour=Depth))  
P <- P + geom_point(size = 8, position=position_dodge(width = 0.75))
P <- P +  geom_errorbar(aes(ymin=proportion_of_variance-proportional_stdev, ymax=proportion_of_variance+proportional_stdev), 
                        width=0.5, size = 2, position=position_dodge(width = 0.75)) + 
  scale_x_discrete(limits=c("between transects", "between plots", "within plots"))
P <- P + theme_bw() + 
  theme(plot.title = element_text(hjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()) + 
  theme(legend.title=element_text()) + labs(x="Source of variance", y="Proportion of soil moisture variance") + 
  theme(legend.position="none")
P <- P + scale_color_manual(values=c("grey", "black"))
P <- P + theme(axis.line = element_line(size = 2))
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
P <- P + scale_y_continuous(limits = c(0, 1), breaks = c(0, .2, .4, .6, .8, 1), labels = scales::percent)
P
