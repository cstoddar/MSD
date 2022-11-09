library(tidyverse)
library(ggplot2)
library(ggpubr)
library(hrbrthemes)
library(viridis)
library(rstatix)

data_full <- read_csv("data_full.csv")

#Convert relevant columns from factors (alpha-numeric categ.) to characters (numeric only categ.)
names <- c('msd_infs', 'msd_moms')
data_full[,names] <- lapply(data_full[,names] , as.character)
str(data_full)

#Create long dataframe with essential data for boxplot
long <- data_full %>%
  rename("229E" = h_co_v_299e_spike,
         "HKU1" = h_co_v_hku1_spike,
         "NL63" = h_cov_nl63_spike,
         "OC43" = h_co_v_oc43_spike,
         "N" = sars_co_v_2_nucleocapsid,
         "RBD" = sars_co_v_2_s1_rbd,
         "NTD" = sars_co_v_2_s1_ntd,
         "Spike" = sars_co_v_2_spike) %>%
  select(msd_testlist, msd_moms, msd_infs, "229E", "HKU1", "NL63", 
         "OC43", "N", "RBD", "NTD", "Spike") %>%
  pivot_longer(
    cols = 4:11,
    names_to = "antigen",
    values_to = "AU_per_mL")
  
#filter for desired sample groups, merge names of sample group columns for moms and infs
long_pre_pan <- long %>%
  filter(msd_moms == "3" | msd_moms == "4" | msd_infs == "3" | msd_infs == "4") %>%
  unite(group_key, msd_infs, msd_moms, na.rm = TRUE)
  
#new boxplot attempt, trying to get p-values on there
bxp2 <- ggboxplot(long_pre_pan, x = "antigen", y = "AU_per_mL",
                  color = "group_key", add = "jitter") +
  scale_x_discrete(
    limits = c("Spike", "RBD", "NTD", "N", "HKU1", "OC43", "229E", "NL63")) #manually order boxes
bxp2
  
#Performing stat test (using rstatix pkg), see https://tinyurl.com/c44x9ttw
stat.test <- long_pre_pan %>%
  group_by(antigen) %>%
  wilcox_test(AU_per_mL ~ group_key) %>%
  adjust_pvalue(method = "fdr") %>% #note it's better to use benjamini-hochberg correction (FDR) 
  add_significance("p.adj")
stat.test

#Add p-values to plot
stat.test <- stat.test %>%
  add_xy_position(x = "antigen", dodge = 0.8) %>%
  mutate(y.position = c(25000, 25000, 25000, 25000, 50000, 112000, 75000, 50000))
bxp2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0) 

#Plotting histograms of groups 3 and 4
hist_filt <- long_pre_pan %>%
  filter(group_key == "4") %>%
  filter(antigen == "OC43") %>%
  mutate(newlog10 = log10(AU_per_mL))

ggplot(data = hist_filt, aes(newlog10)) +
  geom_histogram()


  
