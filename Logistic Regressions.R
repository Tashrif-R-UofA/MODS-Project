pacman::p_load(broom,
               boot)

#Loading the dataset
project_data <- read.delim("MODS_Rohrlach_EAGER.tsv")

## VARIABLE NAME: percent_endogenous_dna_over_30bp ##
project_data <- project_data %>% 
  dplyr:: mutate(percent_endogenous_dna_over_30bp_new = percent_endogenous_dna_over_30bp/100)
#view(project_data$percent_endogenous_dna_over_30bp_new)

# First Glance
ggplot(aes(x = normalized_nr_of_input_reads, y = percent_endogenous_dna_over_30bp_new, col = source, fill = source), data = project_data) +
  theme_bw() +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, linetype = "solid") +
  labs(title = "Scatterplot of proportion of endogenous DNA vs normalized number of input reads for each source",
       x = "Normalized number of input reads",
       y = "Proportion of endogenous DNA")

#Starting with Complex Model
project_data.logr1 <- glm(percent_endogenous_dna_over_30bp_new~normalized_nr_of_input_reads *source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr1)

#Dropping Interaction Term
project_data.logr2 <- glm(percent_endogenous_dna_over_30bp_new~normalized_nr_of_input_reads + source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr2)

# Using WaldTest to find the best model
lmtest::waldtest(project_data.logr1, project_data.logr2)
# Here, p-value 0.7985 is more than 0.05, hence, the models don't explain the data significantly differently, hence, simpler model, additive is the better one.

#Confirming the if additive is the best one
# Making separate models and comparing using WaldTest
project_data.logr1_input <- glm(percent_endogenous_dna_over_30bp_new~normalized_nr_of_input_reads, data = project_data, family = binomial(link = "logit"))
project_data.logr1_source <- glm(percent_endogenous_dna_over_30bp_new~source, data = project_data, family = binomial(link = "logit"))

lmtest::waldtest(project_data.logr2, project_data.logr1_source)
# The p-value is 0.4109>0.05, so, the null hypothesis can't be rejected, that means input_reads can be removed.

lmtest::waldtest(project_data.logr2, project_data.logr1_input)
# The p-value is 0.04107<0.05, so, the null hypothesis can be rejected, that means source cannot be removed

#Visualizing the best model: The model with source
summary(project_data.logr1_source)
glm.diag.plots(project_data.logr1_source)



## VARIABLE NAME: proportion_of_duplicate_reads ##
project_data <- project_data %>% 
  dplyr:: mutate(proportion_of_duplicate_reads_new = proportion_of_duplicate_reads/100)
#view(project_data$proportion_of_duplicate_reads_new)

# First Glance
ggplot(aes(x = normalized_nr_of_input_reads, y = proportion_of_duplicate_reads_new, col = source, fill = source), data = project_data) +
  theme_bw() +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, linetype = "solid") +
  labs(title = "Scatterplot of proportion of duplicate Reads vs normalized number of input reads for each source",
       x = "Normalized number of input reads",
       y = "Proportion of duplicate reads")

#starting with complex model
project_data.logr3 <- glm(proportion_of_duplicate_reads_new~normalized_nr_of_input_reads *source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr3)

#Dropping Interaction Term
project_data.logr4 <- glm(proportion_of_duplicate_reads_new~normalized_nr_of_input_reads + source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr4)

## Using WaldTest to find the best model
lmtest::waldtest(project_data.logr3, project_data.logr4)
# Here, the p-value is more than 0.9919>0.05, hence, the models don't explain the data significantly differently, hence, simpler model, additive is the better one.

# Making separate models and comparing using WaldTest
project_data.logr3_input <- glm(proportion_of_duplicate_reads_new~normalized_nr_of_input_reads, data = project_data, family = binomial(link = "logit"))
project_data.logr3_source <- glm(proportion_of_duplicate_reads_new~source, data = project_data, family = binomial(link = "logit"))

lmtest::waldtest(project_data.logr4, project_data.logr3_source)
# The p-value is 0.973>0.05, so, null hypothesis can't be rejected, that means input_reads can be removed

lmtest::waldtest(project_data.logr4, project_data.logr3_input)
# The p-value is 0.9329>0.05, so, null hypothesis can't, that means source can be removed too

# Comparing source only model with Null Model
project_data.logr3_null <- glm(proportion_of_duplicate_reads_new~1, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr3_null)

lmtest::waldtest(project_data.logr3_source, project_data.logr3_null)
# p-value = 0.9254>0.05, thus null model is the best one.

#Visualizing the best model
summary(project_data.logr3_null)
glm.diag.plots(project_data.logr3_null)


## VARIABLE NAME: percent_gc_of_unique_reads ##
project_data <- project_data %>% 
  dplyr:: mutate(percent_gc_of_unique_reads_new = percent_gc_of_unique_reads/100)
#view(project_data$percent_gc_of_unique_reads_new)

# First Glance
ggplot(aes(x = normalized_nr_of_input_reads, y = percent_gc_of_unique_reads_new, col = source, fill = source), data = project_data) +
  theme_bw() +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, linetype = "solid") +
  labs(title = "Scatterplot of % of covered sites that are Guanine/Cytosine vs normalized number of input reads for each source",
       x = "Normalized number of input reads",
       y = "% of covered sites that are Guanine/Cytosine")

#Starting with Complex Model
project_data.logr5 <- glm(percent_gc_of_unique_reads_new~normalized_nr_of_input_reads *source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr5)

#Dropping Interaction Term
project_data.logr6 <- glm(percent_gc_of_unique_reads_new~normalized_nr_of_input_reads + source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr6)

## Using WaldTest to find the best model
lmtest::waldtest(project_data.logr5, project_data.logr6)
# Here, the p-value is more than 0.981>0.05, hence, the models don't explain the data significantly differently, hence, simpler model, additive is the better one.

# Making separate models and comparing using WaldTest
project_data.logr5_input <- glm(percent_gc_of_unique_reads_new~normalized_nr_of_input_reads, data = project_data, family = binomial(link = "logit"))
project_data.logr5_source <- glm(percent_gc_of_unique_reads_new~source, data = project_data, family = binomial(link = "logit"))

lmtest::waldtest(project_data.logr6, project_data.logr5_source)
# The p-value is 0.9672>0.05, so null hypothesis cannot be rejected, that means input_reads can be removed

lmtest::waldtest(project_data.logr6, project_data.logr5_input)
# The p-value is 0.7781>0.05, so null hypothesis cannot be rejected, that means source can be removed too

# Comparing source only model with Null Model
project_data.logr5_null <- glm(percent_gc_of_unique_reads_new~1, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr5_null)

lmtest::waldtest(project_data.logr5_source, project_data.logr5_null)
# p-value = 0.7797>0.05, thus null model is the best one

#Visualizing the best model
summary(project_data.logr5_null)
glm.diag.plots(project_data.logr5_null)


## VARIABLE NAME: nuclear_contamination_m1_ml ##

# Transform the response variable using absolute value for negative instances
project_data$nuclear_contamination_m1_ml <- ifelse(project_data$nuclear_contamination_m1_ml < 0, 
                                                   abs(project_data$nuclear_contamination_m1_ml), 
                                                   project_data$nuclear_contamination_m1_ml)

# Basic Visualization with transformed response variable
ggplot(aes(x = normalized_nr_of_input_reads, y = nuclear_contamination_m1_ml, col = source, fill = source), data = project_data) +
  theme_bw() +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, linetype = "solid") +
  labs(title = "Scatterplot of logit-transformed measure of nuclear contamination against normalized number of input reads for each source",
       x = "Normalized number of input reads",
       y = "Measure of nuclear contamination")

#Starting with Complex Model
project_data.logr7 <- glm(nuclear_contamination_m1_ml ~ normalized_nr_of_input_reads * source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr7)

#Dropping Interaction Term
project_data.logr8 <- glm(nuclear_contamination_m1_ml~normalized_nr_of_input_reads + source, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr8)

## Using WaldTest to find better model
lmtest::waldtest(project_data.logr7, project_data.logr8)
# Here, the p-value is more than 0.7294>0.05, hence, the models don't explain the data significantly differently, hence, simpler model, additive is the better one.

# Making separate models and comparing using WaldTest
project_data.logr7_input <- glm(nuclear_contamination_m1_ml~normalized_nr_of_input_reads, data = project_data, family = binomial(link = "logit"))
project_data.logr7_source <- glm(nuclear_contamination_m1_ml~source, data = project_data, family = binomial(link = "logit"))

lmtest::waldtest(project_data.logr8, project_data.logr7_source)
# The p-value is 0.8637>0.05, so, null hypothesis cannot be rejected, that means input_reads can be removed

lmtest::waldtest(project_data.logr8, project_data.logr7_input)
# The p-value is 0.8188>0.05, so, null_hypothesis cannot be rejected, that means source can be removed too

# Comparing source only model with Null Model
project_data.logr7_null <- glm(nuclear_contamination_m1_ml~1, data = project_data, family = binomial(link = "logit"))
summary(project_data.logr7_null)

lmtest::waldtest(project_data.logr7_source, project_data.logr7_null)
# p-value = 0.7839>0.05, thus null model is the best one

#Visualizing the best model
summary(project_data.logr7_null)
glm.diag.plots(project_data.logr7_null)