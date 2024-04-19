pacman::p_load(tidyverse, tidymodels, readr, skimr, knitr,
               dplyr, janitor, carplot, vip, forcats, recipes, parsnips, glmnet, tinytex, ggplot2, ggplot)

#Loading the dataset
project_data <- read.delim("MODS_Rohrlach_EAGER.tsv")
--------------------------------------------------------------------------------
##Variable Analysis: mt_to_nuclear_read_ratio

#normalizing the predictor nr_of_input_reads:
mean_nr_of_input_reads <- mean(project_data$nr_of_input_reads)
sd_nr_of_input_reads <- sd(project_data$nr_of_input_reads)
project_data$normalized_nr_of_input_reads <- (
  project_data$nr_of_input_reads - mean_nr_of_input_reads) / sd_nr_of_input_reads

#First Glance
project_data %>%
  ggplot(aes(x=normalized_nr_of_input_reads,
             y=mt_to_nuclear_read_ratio,
             fill=source,
             col=source))+
  theme_bw()+
  geom_point(pch=21)+
  geom_smooth(se=F,
              method='lm')

#Starting with complex linear model
project_data.lm1 <- lm(mt_to_nuclear_read_ratio~normalized_nr_of_input_reads*source,
                       data=project_data)
summary(project_data.lm1)

#Dropping the Interraction Term
project_data.lm2 <- lm(mt_to_nuclear_read_ratio~normalized_nr_of_input_reads+source,
                       data=project_data)
summary(project_data.lm2)

#Dropping non-significant predictor
project_data.lm3 <- lm(mt_to_nuclear_read_ratio~source,
                       data=project_data)
summary(project_data.lm3)

# Comparison using ANOVA
anova(project_data.lm1, project_data.lm3)
# Since p-value = 0.6319>0.05, So, simpler model, project_data.lm3 is the better performing model.

#Visualizing Best Model
par(mfrow = c(2,2))
plot(project_data.lm3)


## Variable Analysis: covered_sn_ps_on_1240k

#normalizing the predictor covered_sn_ps_on_1240k:
mean_covered_sn_ps_on_1240k <- mean(project_data$covered_sn_ps_on_1240k)
sd_covered_sn_ps_on_1240k <- sd(project_data$covered_sn_ps_on_1240k)
project_data$normalized_covered_sn_ps_on_1240k <- (
  project_data$covered_sn_ps_on_1240k - mean_covered_sn_ps_on_1240k) / sd_covered_sn_ps_on_1240k

#First Glance
project_data %>%
  ggplot(aes(x=normalized_nr_of_input_reads,
             y=normalized_covered_sn_ps_on_1240k,
             fill=source,
             col=source))+
  theme_bw()+
  geom_point(pch=21)+
  geom_smooth(se=F,
              method='lm')

#Starting with complex linear model
project_data.lm4 <- lm(normalized_covered_sn_ps_on_1240k~normalized_nr_of_input_reads*source,
                       data=project_data)
summary(project_data.lm4)

# Dropping Interaction Terms
project_data.lm5 <- lm(normalized_covered_sn_ps_on_1240k~normalized_nr_of_input_reads+source,
                       data=project_data)
summary(project_data.lm5)

# Comparison using ANOVA
anova(project_data.lm4, project_data.lm5)
# Since, p-value = 0.07637 >0.05, So, the simpler model, that is, project_data.lm5 is the better performing model.

#Visualizing Best Model
par(mfrow = c(2,2))
plot(project_data.lm5)


## Response Variable: nr_of_mapped_reads_over_30bp

#normalizing the predictor nr_of_mapped_reads_over_30bp:
mean_nr_of_mapped_reads_over_30bp <- mean(project_data$nr_of_mapped_reads_over_30bp)
sd_nr_of_mapped_reads_over_30bp <- sd(project_data$nr_of_mapped_reads_over_30bp)
project_data$normalized_nr_of_mapped_reads_over_30bp <- (
  project_data$nr_of_mapped_reads_over_30bp - mean_nr_of_mapped_reads_over_30bp) / sd_nr_of_mapped_reads_over_30bp

#First Glance
project_data %>%
  ggplot(aes(x=normalized_nr_of_input_reads,
             y=normalized_nr_of_mapped_reads_over_30bp,
             fill=source,
             col=source))+
  theme_bw()+
  geom_point(pch=21)+
  geom_smooth(se=F,
              method='lm')

#Starting with complex linear model
project_data.lm6 <- lm(normalized_nr_of_mapped_reads_over_30bp~normalized_nr_of_input_reads*source,
                       data=project_data)
summary(project_data.lm6)

# Dropping interaction term
project_data.lm7 <- lm(normalized_nr_of_mapped_reads_over_30bp~normalized_nr_of_input_reads+source,
                       data=project_data)
summary(project_data.lm7)

# Comparison using ANOVA
anova(project_data.lm6, project_data.lm7)
# Since, p-value = 0.1636 >0.05, So, the simpler model, that is, project_data.lm7 is the better performing model.

#Visualizing the best model
par(mfrow = c(2,2))
plot(project_data.lm7)


## Variable Analysis: nr_of_unique_mapped_reads

#normalizing the predictor nr_of_unique_mapped_reads:
mean_nr_of_unique_mapped_reads <- mean(project_data$nr_of_unique_mapped_reads)
sd_nr_of_unique_mapped_reads <- sd(project_data$nr_of_unique_mapped_reads)
project_data$normalized_nr_of_unique_mapped_reads <- (
  project_data$nr_of_unique_mapped_reads - mean_nr_of_unique_mapped_reads) / sd_nr_of_unique_mapped_reads

#First Glance
project_data %>%
  ggplot(aes(x=normalized_nr_of_input_reads,
             y=normalized_nr_of_unique_mapped_reads,
             fill=source,
             col=source))+
  theme_bw()+
  geom_point(pch=21)+
  geom_smooth(se=F,
              method='lm')

#Starting with complex linear model
project_data.lm8 <- lm(normalized_nr_of_unique_mapped_reads~normalized_nr_of_input_reads*source,
                       data=project_data)
summary(project_data.lm8)


#Dropping the interaction term
project_data.lm9 <- lm(normalized_nr_of_unique_mapped_reads~normalized_nr_of_input_reads+source,
                       data=project_data)
summary(project_data.lm9)


# Comparison using ANOVA
anova(project_data.lm8, project_data.lm9)
# Since, p-value = 0.2242 >0.05, So, the simpler model, that is, project_data.lm9 is the better performing model.

# Visualizing the best model
par(mfrow = c(2,2))
plot(project_data.lm9)


## Variable Analysis: mean_fold_coverage

#normalizing the predictor mean_fold_coverage:
mean_mean_fold_coverage <- mean(project_data$mean_fold_coverage)
sd_mean_fold_coverage <- sd(project_data$mean_fold_coverage)
project_data$normalized_mean_fold_coverage <- (
  project_data$mean_fold_coverage - mean_mean_fold_coverage) / sd_mean_fold_coverage

#First Glance
project_data %>%
  ggplot(aes(x=normalized_nr_of_input_reads,
             y=normalized_mean_fold_coverage,
             fill=source,
             col=source))+
  theme_bw()+
  geom_point(pch=21)+
  geom_smooth(se=F,
              method='lm')

#Starting with complex linear model
project_data.lm10 <- lm(normalized_mean_fold_coverage~normalized_nr_of_input_reads*source,
                        data=project_data)
summary(project_data.lm10)

#Dropping interaction term
project_data.lm11 <- lm(normalized_mean_fold_coverage~normalized_nr_of_input_reads+source,
                        data=project_data)
summary(project_data.lm11)

# Comparison using ANOVA
anova(project_data.lm10, project_data.lm11)
# Since, p-value = 0.085463  >0.05, So, the simpler model, that is, project_data.lm11 is the better performing model.

#Visualizing the best model
par(mfrow = c(2,2))
plot(project_data.lm11)


## Variable Analysis: nr_sn_ps_used_in_contamination_estimation

#normalizing the predictor nr_sn_ps_used_in_contamination_estimation:
mean_nr_sn_ps_used_in_contamination_estimation <- mean(project_data$nr_sn_ps_used_in_contamination_estimation)
sd_nr_sn_ps_used_in_contamination_estimation <- sd(project_data$nr_sn_ps_used_in_contamination_estimation)
project_data$normalized_nr_sn_ps_used_in_contamination_estimation <- (
  project_data$nr_sn_ps_used_in_contamination_estimation - mean_nr_sn_ps_used_in_contamination_estimation) / sd_nr_sn_ps_used_in_contamination_estimation

#First Glance
project_data %>%
  ggplot(aes(x=normalized_nr_of_input_reads,
             y=normalized_nr_sn_ps_used_in_contamination_estimation,
             fill=source,
             col=source))+
  theme_bw()+
  geom_point(pch=21)+
  geom_smooth(se=F,
              method='lm')

#Starting with Complex model
project_data.lm12 <- lm(normalized_nr_sn_ps_used_in_contamination_estimation~normalized_nr_of_input_reads*source,
                        data=project_data)
summary(project_data.lm12)


# Dropping interaction term
project_data.lm13 <- lm(normalized_nr_sn_ps_used_in_contamination_estimation~normalized_nr_of_input_reads+source,
                        data=project_data)
summary(project_data.lm13)


## Comparison using ANOVA
anova(project_data.lm12, project_data.lm13)
# Since, p-value = 0.003958 < 0.05, So, the simpler model, without the interaction term, project_data.lm13 is the better performing model.

#Visualizing the best model
par(mfrow = c(2,2))
plot(project_data.lm13)

