library(tidyverse)
library(magrittr)
library(funModeling)

source("functions.R")

# input files
inp_dir <- file.path(getwd(), "Input", "fromMarseille")
sample_groups <- read.delim(file.path(inp_dir, "Class_common_all.tsv"), header = FALSE, col.names = c("Sample", "Group"))
peptide <- read.delim(file.path(inp_dir, "PeptidomeMerged_common_all.tsv"), check.names = FALSE)
protein <- read.delim(file.path(inp_dir, "Proteome_common_all.tsv"))
miRNA <- read.delim(file.path(inp_dir, "miRNA_normalized_log2_common_all.tsv"))

# data structure
sum_peptide <- df_status(peptide)
sum_protein <- df_status(protein)
sum_miRNA <- df_status(miRNA)
## colnames are sample_IDs, first col is peptide (internal code), miRNA (hsa-code) or protein ID (uniprot ID)
## peptidome is a 4 decimal number (5 digit), missing are indicated as 0.0000, 7302 in total
## protein is a 4 decimal number (5 digit) scientific number, missing are indicated as 0.0000 e+00, 1613 in total
## miRNA is a 6 decimal number (7 digit), there are no 0 values, 503 in total

# determine current number of samples per group
sample_groups %>% group_by(Group) %>% summarise(n())

# peptide

##data distribution
profiling_num(peptide)
peptide %>%
  pivot_longer(-1) %>%
  filter(value > 0) %>%
  mutate(value = log10(value)) %>%
  # pivot_wider() %>%
  select(-1) %>%
  profiling_num()

peptide %>%
  pivot_longer(-1) %>%
  ggplot(aes(as.factor(.data[[1]]), log10(value), fill = name)) +
    geom_point(pch = 21, alpha = 0.5, stroke = NA) +
    theme_minimal()

peptide %>%
  pivot_longer(-1) %>%
  filter(value > 0) %>%
  ggplot(aes(log10(value))) +
    geom_histogram(bins = 30)

# generate new data
peptide_new <- GenerateDF(peptide, transformation = "log10")
system.time(peptide_new2 <- GenerateDF(peptide, transformation = "log10"))

# assess new data
profiling_num(peptide)
profiling_num(peptide_new)
df_status(peptide)
df_status(peptide_new)

# plot data
peptide_concat <- data.frame(rbind(cbind(peptide %>% pivot_longer(-1), "Source" = "org"),
                                   cbind(peptide_new %>% pivot_longer(-1), "Source" = "new"))) %>%
  left_join(sample_groups, by = c("name" = "Sample")) %>%
  mutate("Group" = ifelse(Source == "new", gsub("_.*", "", name), Group),
         "value" = ifelse(value==0, 0, log10(value)))

peptide_concat %>%
  ggplot(aes(as.factor(.[[1]]), log10(value), colour = Source)) +
  geom_point(shape = 16, alpha = 0.2, stroke = 0, size = 3) +
  theme_minimal()

peptide_concat %>%
  # filter(value > 0) %>%
  ggplot(aes(Group, value)) +
    geom_boxplot(aes(colour = Source), position = position_dodge(0.9)) +
    geom_violin(aes(fill = Source), alpha = 0.4, colour = NA, position = position_dodge(0.9)) +
    theme_minimal()

PlotHist(peptide, peptide_new, transformation = "log10")

df <- peptide %>%
  pivot_longer(-1) %>%
  mutate(value = ifelse(value==0, NA, log10(value))) %>%
  left_join(sample_groups, by = c("name" = "Sample")) %>%
  group_by(.[1], Group) %>%
  summarise("0perc" = sum(is.na(value))/n(),
            "mean" = mean(value, na.rm = TRUE),
            mean = ifelse(is.nan(mean), 0 , mean),
            "sd" = sd(value, na.rm = TRUE),
            sd = ifelse(is.na(sd), 0 , 0.5*sd),
            .groups = "drop")

x <- data.frame(t(apply(df, 1,
  function(x) sample(c(0,10^(rnorm(1, mean = df$mean, sd = df$sd))),
  50, replace = TRUE, prob = c(x[3], 1-x[3])))))
  
  # df2 <- data.frame(t(apply(df, 1, function(x) sample(c(NA, 0), n, replace = TRUE, prob = c(1-x[3], x[3])))))
  
  # for(r in nrow(df)){
  #   
  # }

## Troubleshoot
# set.seed(101)
# df <- data.frame(matrix(sample(c(10:20, 0), 100, replace=TRUE), ncol=10))
# 
# transform(df, "0_perc" = apply(df == 0, 1, sum), "mean" = apply(df, 1, mean), "sd" = apply(df, 1, sd))
#  
# profiling_num(peptide)