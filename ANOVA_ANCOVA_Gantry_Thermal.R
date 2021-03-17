install.packages("psych")
library(psych)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(ggplot2)
library(HH)
plant_data <- read.csv("Desktop/plant_detections_new.csv")
plant_data1 <- read.csv("Desktop/plant_detections_outliers.csv")

#------------------------------------ ANOVAS for Thermal -------------------------------------------
#------------------------------------ Atm Temp -----------------------------------------------------
res.aov2 <- aov(median_c ~ atm_temp, data = plant_data)
summary(res.aov2)
plot(res.aov2)
describeBy(plant_data$median_c, plant_data$atm_temp)
#------------------------------------ Treatment ----------------------------------------------------
res.aov2 <- aov(median_c ~ treatment, data = plant_data)
summary(res.aov2)
plot(res.aov2)
describeBy(plant_data$median_c, plant_data$treatment)
#------------------------------------ Wind velocity ------------------------------------------------
res.aov2 <- aov(median_c ~ azmet_wind_velocity, data = plant_data)
summary(res.aov2)
plot(res.aov2)
describeBy(plant_data$median_c, plant_data$azmet_wind_velocity)
#------------------------------------ Genotype -----------------------------------------------------
res.aov2 <- aov(median_c ~ genotype, data = plant_data)
summary(res.aov2)
describeBy(plant_data$median_c, plant_data$genotype)
#--------------------------------- ANCOVAS for thermal ---------------------------------------------
#---------------------------------------------------------------------------
ancova<-lm(median_c~atm_temp*treatment, data=plant_data)
anova(ancova)
summary(ancova)
plot(ancova)

t1<-plant_data[which(plant_data$treatment=='treatment 1'),]
t2<-plant_data[which(plant_data$treatment=='treatment 2'),]
t3<-plant_data[which(plant_data$treatment=='treatment 3'),]
border <- plant_data[which(plant_data$treatment=='border'),]
#----------------------------- Treatment 1 ----------------------------------------------  
attach(t1)
plot(t1$atm_temp, t1$median_c, xlab="atm_temp", ylab="PCT")
abline(lm(median_c~atm_temp), col="red")
fit<-lm(median_c~atm_temp, data=t1)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)
#----------------------------- Treatment 2 ----------------------------------------------  
attach(t2)
plot(t2$atm_temp, t2$median_c, xlab="atm_temp", ylab="PCT", pch=19)
abline(lm(median_c~atm_temp), col="red")
fit<-lm(median_c~atm_temp, data=t2)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)
#----------------------------- Treatment 3 ----------------------------------------------  
attach(t3)
plot(t3$atm_temp, t3$median_c, xlab="atm_temp", ylab="PCT", pch=19)
abline(lm(median_c~atm_temp), col="red")
fit<-lm(median_c~atm_temp, data=t3)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)
#----------------------------- Border ---------------------------------------------- 
attach(border)
plot(border$atm_temp, border$median_c, xlab="atm_temp", ylab="PCT", pch=19)
abline(lm(median_c~atm_temp), col="red")
fit<-lm(median_c~atm_temp, data=border)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)

#----------------------------- Presenting Results with Azmet Temp----------------------------------------
ggplot(plant_data, (aes(x=atm_temp, y=median_c, color=treatment, shape=treatment))) + theme_classic() +ylab("PCT") + 
  xlab("atm_temp") + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

ggplot(plant_data, (aes(x=atm_temp, y=median_adj))) + theme_classic() +ylab("Adj PCT") + 
  xlab("atm_temp") + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)
#----------------------------- Running ANCOVA using Env Temperature ----------------------------------------

ancova<-lm(median_c~Env_Temperature*treatment, data=plant_data)
anova(ancova)
summary(ancova)

t1<-plant_data[which(plant_data$treatment=='treatment 1'),]
t2<-plant_data[which(plant_data$treatment=='treatment 2'),]
t3<-plant_data[which(plant_data$treatment=='treatment 3'),]
border <- plant_data[which(plant_data$treatment=='border'),]
#----------------------------- Treatment 1 ----------------------------------------------  
attach(t1)
plot(t1$Env_Temperature, t1$median_c, xlab="env_temp", ylab="PCT")
abline(lm(median_c~Env_Temperature), col="red")
fit<-lm(median_c~Env_Temperature, data=t1)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)
#----------------------------- Treatment 2 ----------------------------------------------  
attach(t2)
plot(t2$Env_Temperature, t2$median_c, xlab="env_temp", ylab="PCT")
abline(lm(median_c~Env_Temperature), col="red")
fit<-lm(median_c~Env_Temperature, data=t2)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)
#----------------------------- Treatment 3 ----------------------------------------------  
attach(t3)
plot(t3$Env_Temperature, t3$median_c, xlab="env_temp", ylab="PCT")
abline(lm(median_c~Env_Temperature), col="red")
fit<-lm(median_c ~ Env_Temperature, data=t3)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)
#----------------------------- Border ---------------------------------------------- 
attach(border)
plot(border$Env_Temperature, border$median_c, xlab="env_temp", ylab="PCT")
abline(lm(median_c~Env_Temperature), col="red")
fit<-lm(median_c~Env_Temperature, data=border)
library(stats)
plot(density(fit$residuals))
qqnorm(fit$residuals)
qqline(fit$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))
summary(fit)

#----------------------------- Presenting Results with Env Temp ----------------------------------------
ggplot(plant_data, (aes(x=Env_Temperature, y=median_c, color=treatment, shape=treatment))) + theme_classic() +ylab("PCT") + 
  xlab("atm_temp") + geom_boxplot()

ggplot(plant_data, (aes(x=atm_temp, y=median_adj2, color=treatment, shape=treatment))) + theme_classic() +ylab("Adj PCT") + 
  xlab("atm_temp") + geom_boxplot()

#--------------------------------Two way ANCOVA--------------------------------------
res.aov <- plant_data %>% anova_test(median_c ~ atm_temp + treatment*relative_humidity)
get_anova_table(res.aov)


# Effect of treatment at each level of exercise
plant_data %>%
  group_by(treatment) %>%
  anova_test(median_c ~ atm_temp + relative_humidity + azmet_wind_velocity)

# Pairwise comparisons
library(emmeans)
pwc <- plant_data %>%
  group_by(treatment) %>%
  emmeans_test(
    median_c ~ treatment, covariate = atm_temp, 
    p.adjust.method = "bonferroni"
  )
pwc 

# Shows the estimated means
get_emmeans(pwc)

# Line plot
ggline(
  get_emmeans(pwc), x = "treatment", y = "emmean", color = "treatment", palette = "jco") +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high, color = treatment), 
    width = 0.1
  )
# Normal means calculated from the data
library(plyr)
normal_mean = ddply(plant_data, .(treatment), summarize,  median_mean=mean(median_c))
df1 = data.frame(normal_mean)

ggplot(data=df1, aes(x=treatment, y=median_mean, group=1)) + geom_line(color = "red") + geom_point() + 
  geom_line(data = get_emmeans(pwc), aes(x=treatment, y = emmean), color = "blue") +
  geom_errorbar(data = get_emmeans(pwc), aes(ymin = conf.low, ymax = conf.high, color = treatment), width = 0.1)

#---------------------------------- One Way ANCOVA --------------------------------------

res.aov <- plant_data %>% anova_test(median_c ~ atm_temp + treatment)
get_anova_table(res.aov)
plot(res.aov)

# Pairwise comparisons
library(emmeans)
pwc <- plant_data %>% 
  emmeans_test(
    median_c ~ treatment, covariate = atm_temp,
    p.adjust.method = "bonferroni"
  )
pwc

# Shows the estimated means
get_emmeans(pwc)

# Plotting new estimated means by treatment
pwc <- pwc %>% add_xy_position(x = "treatment", fun = "mean_se")
ggline(get_emmeans(pwc), x = "treatment", y = "emmean") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
  stat_pvalue_manual(pwc, hide.ns = TRUE, tip.length = FALSE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

library(plyr)
# Calculate the mean of the treatments before adjusting the mean
normal_mean = ddply(plant_data, .(treatment), summarize,  median_mean=mean(median_c))
df = data.frame(normal_mean)
# Graph both the adjusted and the normal mean together
ggplot(data=df, aes(x=treatment, y=median_mean, group=1)) +
  geom_line(color = "red") + geom_point() + geom_line(data = get_emmeans(pwc), aes(x=treatment, y = emmean), color = "blue") + 
  ggtitle ("2020-03-03") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + 
  xlab("Irrigation Treatment")

#-----------------------------------Using ANCOVA for adjusting the means------------------------------------------
model.1 = lm (median_c ~ atm_temp + treatment + relative_humidity + azmet_wind_velocity, data = plant_data)
Anova(model.1, type="III")
plot(model.1)

model.2 = lm (median_c ~ atm_temp + treatment + relative_humidity, data = plant_data)
Anova(model.1, type="III")

model.3 = lm (median_c ~ atm_temp + treatment, data = plant_data)
Anova(model.1, type="III")

library(effects)
# se=TRUE: show standard errors
adjustedMeans<-effect("treatment", model.1, se=TRUE)
adj_df1 <- data.frame(adjustedMeans)
names(adj_df1)[names(adj_df1) == 'fit'] <- 'adj1'
adj_df1

adjustedMeans<-effect("treatment", model.2, se=TRUE)
adj_df2 <- data.frame(adjustedMeans)
names(adj_df2)[names(adj_df2) == 'fit'] <- 'adj2'
adj_df2

adjustedMeans<-effect("treatment", model.3, se=TRUE)
adj_df3 <- data.frame(adjustedMeans)
names(adj_df3)[names(adj_df3) == 'fit'] <- 'adj3'
adj_df3

library(plyr)
normal_mean = ddply(plant_data, .(treatment), summarize,  median_mean=mean(median_c))
norm_df = data.frame(normal_mean)
norm_df

merged_df1 <- merge(x = adj_df1, y = norm_df, by = "treatment", all = TRUE)

merged_df2 <- merge(x = adj_df2, y = adj_df3, by = "treatment", all = TRUE)

merged_df3 <- merge(x = merged_df1, y = merged_df2, by = "treatment", all = TRUE)
merged_df3

ggplot(data = merged_df3, aes(x = treatment, y = median_mean, group = 1)) + geom_line(color = "black") + 
  geom_line(data = merged_df3, aes(x = treatment, y = adj1, group = 1), color = "blue") + geom_point() + 
  geom_line(data = merged_df3, aes(x = treatment, y = adj2, group = 1), color = "green") + 
  geom_line(data = merged_df3, aes(x = treatment, y = adj3, group = 1), color = "red") 


#-----------------------Using ANCOVA for adjusting the means excluding influential points-------------------------
model.1 = lm (median_c ~ atm_temp + treatment + relative_humidity + azmet_wind_velocity, data = plant_data1)
Anova(model.1, type="III")
plot(model.1)

model.2 = lm (median_c ~ atm_temp + treatment + relative_humidity, data = plant_data1)
Anova(model.2, type="III")

model.3 = lm (median_c ~ atm_temp + treatment, data = plant_data1)
Anova(model.3, type="III")

library(effects)
# se=TRUE: show standard errors
adjustedMeans<-effect("treatment", model.1, se=TRUE)
adj_df1 <- data.frame(adjustedMeans)
names(adj_df1)[names(adj_df1) == 'fit'] <- 'adj1'
adj_df1

adjustedMeans<-effect("treatment", model.2, se=TRUE)
adj_df2 <- data.frame(adjustedMeans)
names(adj_df2)[names(adj_df2) == 'fit'] <- 'adj2'
adj_df2

adjustedMeans<-effect("treatment", model.3, se=TRUE)
adj_df3 <- data.frame(adjustedMeans)
names(adj_df3)[names(adj_df3) == 'fit'] <- 'adj3'
adj_df3

library(plyr)
normal_mean = ddply(plant_data1, .(treatment), summarize,  median_mean=mean(median_c))
norm_df = data.frame(normal_mean)
norm_df

merged_df1 <- merge(x = adj_df1, y = norm_df, by = "treatment", all = TRUE)

merged_df2 <- merge(x = adj_df2, y = adj_df3, by = "treatment", all = TRUE)

merged_df3 <- merge(x = merged_df1, y = merged_df2, by = "treatment", all = TRUE)
merged_df3
#write.csv(merged_df3, 'Adjusted_Means.csv')


ggplot(data = merged_df3, aes(x = treatment, y = median_mean, group = 1)) + geom_line(color = "black") + 
  geom_line(data = merged_df3, aes(x = treatment, y = adj1, group = 1), color = "blue") + 
  geom_line(data = merged_df3, aes(x = treatment, y = adj2, group = 1), color = "green") + 
  geom_line(data = merged_df3, aes(x = treatment, y = adj3, group = 1), color = "red") + 
  ylab('Median Plant Canopy Temperature (°C)') + xlab('Irrigation Treatment') + 
  ggtitle('Adjusted PCT Means by Treatment') + theme(plot.title = element_text(hjust = 0.5))


#---------------------------------------------------------------------
model.1 = lm (median_c ~ atm_temp + genotype + treatment + relative_humidity + azmet_wind_velocity, data = plant_data)
Anova(model.1, type="III")


library(effects)
# se=TRUE: show standard errors
adjustedMeans<-effect("treatment", model.1, se=TRUE)
adj_df1 <- data.frame(adjustedMeans)
names(adj_df1)[names(adj_df1) == 'fit'] <- 'adj1'
adj_df1

library(plyr)
normal_mean = ddply(plant_data1, .(treatment), summarize,  median_mean=mean(median_c))
norm_df = data.frame(normal_mean)
norm_df

merged_df1 <- merge(x = adj_df1, y = norm_df, by = "treatment", all = TRUE)
merged_df1

ggplot(data = merged_df1, aes(x = treatment, y = median_mean, group = 1)) + geom_line(color = "black") + 
  geom_line(data = merged_df1, aes(x = treatment, y = adj1, group = 1), color = "blue")




