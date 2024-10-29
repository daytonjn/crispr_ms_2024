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
#figure 2 plot with effect of mutation on rhythmic strength
#effect of mutation on rhythmic strength (difference in proportion)
rs_df <- read_csv("./data/rhythmic_strength_prop.csv",
                  col_types = cols(treatment =
                                     col_factor(levels = c("control", "period", "pdfr")),
                                   injection_period = col_factor(levels = c("2022", "2023"))))

rs_df$diff[which(rs_df$diff==0)] = 0.005
p_rs = ggplot(data = rs_df,
              aes(x = treatment, y = diff, fill = treatment)) +
  geom_col(position = position_dodge(), alpha = 0.5) +
  geom_errorbar(aes(ymin = diff_lb, ymax = diff_ub), width=0.2, position = position_dodge(0.9))+
  facet_wrap(~injection_period)+
  theme_bw() +
  labs(x = "", y = "Difference in Rhythmic Strength") +
  scale_fill_manual(values = c('period' = "#F8766D",
                               'control' = "#00BFC4",
                               'pdfr' = 'darkorchid'),
                    labels = c('control' = "Wildtype", 
                               'pdfr' = expression(italic(pdfr)~F[0]),
                               'period' = expression(italic(period)~F[0]))) +
  ylim(-0.5,0.4)  +
  scale_x_discrete(labels = c(expression(Wildtype),
                              expression(italic(period)~F[0]),
                              expression(italic(pdfr)~F[0])), name = "") +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.5) +
  geom_text(data = m1_multcomp_b,
            aes(label = .group, x = treatment, y = rep(0.35, 6)),
            color = 'black', size = 4) + facet_wrap(~injection_period)

p_rs = ggplot(data = rs_df,
              aes(x = treatment, y = diff, fill = treatment)) +
  geom_col(position = position_dodge(), alpha = 0.5) +
  geom_errorbar(aes(ymin = diff_lb, ymax = diff_ub), width=0.2, position = position_dodge(0.9))+
  facet_wrap(~injection_period)+
  theme_bw() +
  labs(x = "", y = "Difference in Rhythmic Strength") +
  scale_fill_manual(values = c('period' = "#F8766D",
                               'control' = "#00BFC4",
                               'pdfr' = 'darkorchid'),
                    labels = c('control' = "Wildtype", 
                               'pdfr' = expression(italic(pdfr)~F[0]),
                               'period' = expression(italic(period)~F[0]))) +
  ylim(-0.5,0.4)  +
  scale_x_discrete(labels = c(expression(Wildtype),
                              expression(italic(period)~F[0]),
                              expression(italic(pdfr)~F[0])), name = "") +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 0.5) +
  geom_text(data = m1_multcomp_b,
            aes(label = .group, x = treatment, y = rep(0.35, 6)),
            color = 'black', size = 4) + facet_wrap(~injection_period)

ggsave("./plots/figure2.png", plot = p_rs, width = 8, height = 4, dpi = 800, bg = 'white')

