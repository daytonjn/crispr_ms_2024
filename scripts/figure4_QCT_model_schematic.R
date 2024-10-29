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


######################
#figure 4 model

#plot for modified QCT approach
#design1 (A1W AND A2W)
design1 = data.frame(allele_x = factor(c('1', '1.5', '1', '1.5')),
                     treatment = factor(c("control", "control", 'mutant', 'mutant')),
                     phenotype_y = c(1,2,2,10))

#design2 (A1A2+A1W AND A1A2+A2W)
design2 = data.frame(allele_x = factor(c('1', '2', '1', '2')),
                     treatment = factor(c("control", "control", 'mutant', 'mutant')),
                     phenotype_y = c(1.25,1.75,4,8))



#design1 is with known genotypes
p_d1 = ggplot(data=design1 %>% mutate(treatment = factor(treatment, levels = c('mutant', 'control')))) + 
  geom_line(aes(x = allele_x, y = phenotype_y, colour = treatment, group=treatment), linewidth = 1.5) +
  scale_colour_manual(values = c('mutant' = "orange",
                                 'control' = "#00BFC4"),
                      labels = c('mutant' = "Mutant (Target sgRNA)",
                                 'control' = "Wildtype (Control sgRNA)")) +
  scale_y_continuous(name = "Phenotype (Trait value)",
                     limits = c(0,12)) +
  scale_x_discrete(name = " ",
                   labels = c(expression(A[1]~W), expression(A[2]~W)),
                   expand = c(0.1, 0.1)) +
  theme_bw()  +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25,0.80),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p_d1

#design2 is with unknown genotypes/heterozygotes included
p_d2 = ggplot(data=design2 %>% mutate(treatment = factor(treatment, levels = c('mutant', 'control')))) + 
  geom_line(aes(x = allele_x, y = phenotype_y, colour = treatment, group=treatment), linewidth = 1.5) +
  scale_colour_manual(values = c('mutant' = "orange",
                                 'control' = "#00BFC4"),
                      labels = c('mutant' = "Mutant (Target sgRNA)",
                                 'control' = "Wildtype (Control sgRNA)")) +
  scale_y_continuous(name = "Phenotype (Trait value)",
                     limits = c(0,12)) +
  scale_x_discrete(name = "Genotype (Candidate locus)",
                   labels = c(expression(A[1]~A[2]~+A[1]~W), expression(A[1]~A[2]~+A[2]~W)),
                   expand = c(0.1, 0.1)) +
  theme_bw()  +
  theme(text = element_text(size = 14, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25,0.80),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())


p_d2
p_model = ggarrange(p_d1, p_d2, nrow=2, common.legend = TRUE)

p_model


#strong evidence with interaction by genotype vs. no interaction by genotype

#plot for modified QCT approach
#additive effect with no background interaction
design1 = data.frame(allele_x = factor(c('1', '1.5', '1', '1.5')),
                     treatment = factor(c("control", "control", 'mutant', 'mutant')),
                     phenotype_y = c(1.5,2.5,5.5,6.5))

#additive effect with a background interaction
design2 = data.frame(allele_x = factor(c('1', '2', '1', '2')),
                     treatment = factor(c("control", "control", 'mutant', 'mutant')),
                     phenotype_y = c(1.5,2.5,4.75,10.25))



#design1 is with known genotypes
p_d1 = ggplot(data=design1 %>% mutate(treatment = factor(treatment, levels = c('mutant', 'control')))) + 
  geom_line(aes(x = allele_x, y = phenotype_y, colour = treatment, group=treatment), linewidth = 1.5) +
  scale_colour_manual(values = c('mutant' = "orange",
                                 'control' = "#00BFC4"),
                      labels = c('mutant' = "Mutant (Target sgRNA)",
                                 'control' = "Wildtype (Control sgRNA)")) +
  scale_y_continuous(name = "Phenotype (Trait value)",
                     limits = c(0,12)) +
  scale_x_discrete(name = "Genotype (Candidate locus)",
                   labels = c(expression(A[1]~A[2]~+A[1]~W), expression(A[1]~A[2]~+A[2]~W)),
                   expand = c(0.1, 0.1)) +
  theme_bw()  +
  theme(text = element_text(size = 12, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25,0.80),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p_d1

#design2 is with unknown genotypes/heterozygotes included
p_d2 = ggplot(data=design2 %>% mutate(treatment = factor(treatment, levels = c('mutant', 'control')))) + 
  geom_line(aes(x = allele_x, y = phenotype_y, colour = treatment, group=treatment), linewidth = 1.5) +
  scale_colour_manual(values = c('mutant' = "orange",
                                 'control' = "#00BFC4"),
                      labels = c('mutant' = "Mutant (Target sgRNA)",
                                 'control' = "Wildtype (Control sgRNA)")) +
  scale_y_continuous(name = "Phenotype (Trait value)",
                     limits = c(0,12)) +
  scale_x_discrete(name = " ",
                   labels = c(expression(A[1]~A[2]~+A[1]~W), expression(A[1]~A[2]~+A[2]~W)),
                   expand = c(0.1, 0.1)) +
  theme_bw()  +
  theme(text = element_text(size = 12, family = "sans", color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25,0.80),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())


p_d2
p_model = ggarrange(p_d2, p_d1, nrow=2, common.legend = TRUE)

p_model


library('png')
img1 <- readPNG("./p_model.png")
im_A <- ggplot() +
  background_image(img1) +
  theme(plot.background = element_rect(color = 'white'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t=0.25, l=0.7, r=0.3, b=0.5, unit = "in"))

im_B = ggplot() +
  theme(plot.background = element_rect(color = 'white'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t=0.25, l=0, r=0.25, b=0.5, unit = "in"))

model_plots = ggarrange(p_d2, p_d1, nrow=2, common.legend = TRUE, legend = "top") %>%
  annotate_figure(left = text_grob("Phenotype (Trait value)", color='black', rot=90, size=14))


temp = ggarrange(im_B, model_plots, ncol=2, labels = c("A", "B"), widths = c(1,0.75))
temp
ggsave("./plots/figure4.png", plot = temp, width = 10, height = 4, dpi = 500, units=c('in'),
       bg='white')
