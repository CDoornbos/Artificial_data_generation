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

## data distribution
profiling_num(peptide)

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

## new data
set.seed(101)
peptide_log10 <- peptide %>%
  column_to_rownames(colnames(.)[1]) %>%
  mutate(across(starts_with("C"), function(x) ifelse(x == 0, NA, log10(x)))) %>%
  mutate("0_perc" = apply(is.na(.), 1, sum)/ncol(.),
         "mean" = apply(., 1, function(x) mean(x, na.rm = TRUE)),
         "mean" = ifelse(is.nan(mean), 0, mean),
         "sd" = apply(., 1, function(x) sd(x, na.rm = TRUE)),
         "sd" = ifelse(is.na(sd), 0, sd),
         .keep = "used")

peptide_art <- data.frame(t(apply(peptide_log10, 1, function(x) rnorm(50, mean = x[2], sd = x[3])))) %>%
  mutate(across(is.numeric, function(x){
    10^x
    format(x, scientific = FALSE)
    round(x, digits = 4)
  }))

## missing data distribution
sum_peptide[2:32,] %>%
  ggplot(aes(p_zeros)) +
    geom_histogram(bins = 7)

## add missing data points # TODO here: fix distributions
set.seed(101)
peptide_art <- data.frame(t(apply(peptide_log10, 1, function(x) sample(c(NA, 0), 50, replace = TRUE, prob = c(1-x[1], x[1]))))) %>%
  mutate(across(is.numeric, function(x) ifelse(is.na(x), 10^(rnorm(1, mean = peptide_log10$mean, sd = peptide_log10$sd)), x)))

profiling_num(peptide_art)
df_status(peptide_art)





peptide_new <- GenerateDF(peptide, transformation = "log10")
peptide_new <- format(peptide_new, scientific = FALSE) %>%
  mutate(across(where(is.character), as.numeric))

    
profiling_num(peptide)
profiling_num(peptide_new)
df_status(peptide)
df_status(peptide_new)

peptide_concat <- rbind(cbind(peptide %>% pivot_longer(-1), "Source" = "org"),
                        cbind(peptide_new %>% pivot_longer(-1), "Source" = "new"))

peptide_concat %>%
  ggplot(aes(as.factor(.[[1]]), log10(value), colour = Source)) +
  geom_point(shape = 16, alpha = 0.2, stroke = 0, size = 3) +
  theme_minimal()

peptide_concat %>%
  filter(value > 0) %>%
  ggplot(aes(log10(value), fill = Source)) +
  geom_histogram(bins = 30, alpha = 0.5)


rbind(cbind(peptide %>% pivot_longer(-1), "Source" = "org"),
      cbind(peptide_new %>% pivot_longer(-1), "Source" = "new")) %>%
  filter(value > 0) %>%
  mutate(value = log10(value)) %>%
  group_by(cut(value, 30), Source) %>%
  summarise("number" = n(), .groups = "drop") %>%
  mutate("value" = gsub(".(.*)]", "\\1", .[[1]])) %>%
  separate(value, into = c("from", "to"), sep = ",", convert = TRUE) %>%
  rowwise() %>%
  mutate("mean" = mean(c(from, to), na.rm = TRUE)) %>%
  mutate("number" = ifelse(Source == "new", number*(ncol(peptide)/ncol(peptide_new)), number)) %>%
  ggplot(aes(mean, number, fill = Source)) +
    geom_col(position = "dodge") +
    theme_minimal()


# peptide_art <- data.frame(t(apply(peptide_log10, 1, function(x) sample(c(NA, 0), 50, replace = TRUE, prob = c(1-x[1], x[1]))))) %>%
#   mutate(across(is.numeric, function(x) ifelse(is.na(x), 10^(rnorm(1, mean = peptide_log10$mean, sd = peptide_log10$sd)), x)))
  
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