library(reshape2)
library(readr)
library(tidyverse)
library(bbmle)      # for ICtab()
library(pscl)       # for zi models
library(car)        # for Anova()
library(lme4)       # for lmer() and glmer() 
library(emmeans)    # for multiple comparisons!
library(multcomp)
library(multcompView)
library(ggpubr)
library(lmtest)
library(report)

###########################
#mating success
treatment <- c("period1","period2","pdfr","period1","period2","pdfr")
year <- c(2022,2022,2022,2023,2023,2023)
success <- c(25,20,7,17,22,9)
fail <- c(26,22,2,10,21,5)
ID <- c(0,1,2,0,1,2)
df1 <- data.frame(treatment,year,success, fail,ID)

df1$treatment = factor(df1$treatment)
df1$year = factor(df1$year)

m1_mating <- glm(cbind(success,fail) ~ treatment+year, 
              data = df1, 
              family = binomial(link = "identity"))


car::Anova(m1_mating, type = "2")

#somatic editing
treatment <- c("period1","period2","pdfr","period1","period2","pdfr")
year <- c(2022,2022,2022,2023,2023,2023)
success <- c(11,12,5,14,16,3)
fail <- c(8,3,2,1,3,2)
ID <- c(0,1,2,0,1,2)
df2 <- data.frame(treatment,year,success, fail,ID)

df2$treatment = factor(df2$treatment)
df2$year = factor(df2$year)

m2_somatic <- glm(cbind(success,fail) ~ treatment*year, 
                 data = df2, 
                 family = binomial(link = "identity"))

summary(m2_somatic)
car::Anova(m2_somatic, type = "2")


df3 <- read_csv("./data/indel_percentage.csv", 
                             col_types = cols(injection_period = col_factor(levels = c("2022","2023")),
                                              treatment = col_factor(levels = c("period3", "period4", "pdfr")), indel_percentage = col_number()))

df3$indel_percentage = df3$indel_percentage*100
wilcox.test(df3$indel_percentage[df3$injection_period==2022], df3$indel_percentage[df3$injection_period==2023])

m3_indel <- lm(indel_percentage ~ injection_period*treatment, df3)
pairs(emmeans(m3_indel, ~ injection_period | treatment,type = "response"))

emmeans(m3_indel, ~ injection_period,type = "response")

t.test(indel_percentage ~ injection_period, data = df3)

t.test(indel_percentage ~ injection_period, data = df3[df3$treatment!='pdfr',])


#indel characteristics
df4 <- read_csv("./data/indel_mutations.csv", 
                            col_types = cols(injection_period = col_factor(levels = c("2022", "2023")), 
                                             treatment = col_factor(levels = c("period_sgrna1","period_sgrna2", "pdfr_sgrna"))))

m4_indel <- lm(indel ~ treatment, df4 %>% filter(indel<0))
Anova(m4_indel)
emmeans(m4_indel, ~ treatment,type = "response")

df4_res = df4 %>% group_by(treatment) %>%
  summarize(n = n(),
            mu = mean(indel),
            mode = mode(indel),
            del_n = sum(indel<0),
            ins_n = sum(indel>0),
            del_percent = sum(indel<0)/n,
            inframe_n = sum(indel %% 3 == 0),
            frameshift_n = sum(indel %% 3 != 0),
            frameshift_percent = sum(indel %% 3 != 0)/n)
df4_res
t.test(c(0.818,0.64,1))#deletions
t.test(c(0.818,0.8,0.714)) #frameshift

t.test(c(0.81,0.65,0.52)) #survival to eclosion
t.test(c(0.49,0.48,0.78,0.63,0.51,0.64)) #mating success
t.test(c(0.58,0.80,0.71,0.93,0.84,0.6)) #somatic edits
t.test(c(0.86,1,1,1,1,1)) #somatic edits

#average live adults/cluster of eggs
#67, 50, 11
#12, 12, 4

#56, 79, 33
#8, 8, 4
t.test(c(67/12, 50/12, 11/4, 56/8, 79/8, 33/4))

c(67/12, 50/12, 11/4, 56/8, 79/8, 33/4)
t.test(c(5.58,4.2,2.75),c(7,9.87,8.25))

del_n <- c(16,9)
res <- chisq.test(del_n, p = c(1/2,1/2))
res

df5_res = df4  %>%
  summarize(n = n(),
            mu = mean(indel),
            mode = mode(indel),
            del_n = sum(indel<0),
            ins_n = sum(indel>0),
            del_percent = sum(indel<0)/n,
            inframe_n = sum(indel %% 3 == 0),
            frameshift_n = sum(indel %% 3 != 0),
            frameshift_percent = sum(indel %% 3 != 0)/n)

chisq.test(c(41,13), p = c(1/2,1/2)) #expected 50:50 del:ins
chisq.test(c(43,11), p = c(2/3,1/3)) #expected 66:33 frameshift:inframe

#2022
chisq.test(c(18,4), p = c(2/3,1/3)) #expected 66:33 frameshift:inframe
#2023
chisq.test(c(23,9), p = c(2/3,1/3)) #expected 66:33 frameshift:inframe

#calculate
#proportion of somatic mutants
61/80
n = 80
# Get the proportion
p_hat = 61/n
alpha=0.05

# Calculate the critical z-score
z = qnorm(1-alpha/2)

# Compute the CI
p_hat + c(-1,1)*z*sqrt(p_hat*(1-p_hat)/n)




