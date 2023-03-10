library(tidyverse)
library(magrittr)
library(pheatmap)

set.seed(101)

GenerateDF <- function(df, groups = sample_groups, n = 25, transformation = FALSE, scientific = FALSE, digits = 4){
  if(transformation == FALSE){
    df %<>%
      pivot_longer(-1)
  } else if(transformation == "log10"){
    df %<>%
      pivot_longer(-1) %>%
      mutate(value = ifelse(value==0, NA, log10(value)))
  } else{
    warningCondition("The value of 'transformation' should be FALSE or log10")
  }
  
  #get stats
  df %<>%
    left_join(groups, by = c("name" = "Sample")) %>%
    group_by(.[1], Group) %>%
    summarise("0perc" = sum(is.na(value))/n(),
              "mean" = mean(value, na.rm = TRUE),
              mean = ifelse(is.nan(mean), 0 , mean),
              "sd" = sd(value, na.rm = TRUE),
              sd = ifelse(is.na(sd), 0 , 0.5*sd),
              .groups = "drop")
  
  #make df missing data
  df_new <- c()
  for (group in levels(as.factor(groups$Group))){
    df_split <- df %>%
      filter(Group == group) %>%
      select(-Group)
    df_grp <- apply(df_split, 1, function(x) sample(c(NA, 0), n, replace = TRUE, prob = c(1-x[2], x[2])))
    rownames(df_grp) <- sapply(1:25, function(x) paste0(group, "_", x))
    df_new <- rbind(df_new, cbind(df_grp))
  }
  
  df_new <- data.frame(unique(df[1]), t(df_new), check.names = FALSE)
  
  #impute data
  for (r in 1:nrow(df_new)){
    for (c in 1:ncol(df_new)){
      if (is.na(df_new[r,c])){
        df_sub <- df %>%
          filter(Group == gsub("_.*", "", colnames(df_new[c]))) %>%
          filter(.[[1]] == df_new[r,1])
        if(transformation == FALSE){
          df_new[r,c] <- round(rnorm(1, mean = df_sub$mean, sd = df_sub$sd), digits = digits)
        } else if(transformation == "log10"){
          df_new[r,c] <- round(10^(rnorm(1, mean = df_sub$mean, sd = df_sub$sd)), digits = digits)
        }
      }
    }
  }
  
  #correct format
  df_new <- format(df_new, scientific = scientific)
  
  return(df_new)
}




## alternative function that might be quicker
GenerateDF <- function(df, groups = sample_groups, n = 25, transformation = FALSE, scientific = FALSE, digits = 4){
  #remove
  df = peptide
  groups = sample_groups
  n = 25
  transformation = "log10"
  scientific = FALSE
  digits = 4
  
  if(transformation == FALSE){
    df %<>%
      pivot_longer(-1)
  } else if(transformation == "log10"){
    df %<>%
      pivot_longer(-1) %>%
      mutate(value = ifelse(value==0, NA, log10(value)))
  } else{
    warningCondition("The value of 'transformation' should be FALSE or log10")
  }
  
  #get stats
  df %<>%
    left_join(groups, by = c("name" = "Sample")) %>%
    group_by(.[1], Group) %>%
    summarise("0perc" = sum(is.na(value))/n(),
              "mean" = mean(value, na.rm = TRUE),
              mean = ifelse(is.nan(mean), 0 , mean),
              "sd" = sd(value, na.rm = TRUE),
              sd = ifelse(is.na(sd), 0 , 0.5*sd),
              .groups = "drop")
  
  #make df #TODO: here "Error in sample.int(length(x), size, replace, prob):incorrect number of probabilities"
  
  df_new <- c()
  for (group in levels(as.factor(groups$Group))){
    df_split <- df %>%
      filter(Group == group) %>%
      select(-Group)
    df_grp <- apply(df_split, 1,
                    function(x) sample(c(0,
                      sample(c(0,10^(rnorm(1, mean = x[3], sd = x[4]))))),
                      # NA),
                      n, replace = TRUE, prob = c(x[2], 1-x[2])))
    rownames(df_grp) <- sapply(1:25, function(x) paste0(group, "_", x))
    df_new <- rbind(df_new, cbind(df_grp))
  }
  
  df_new <- data.frame(unique(df[1]), t(df_new), check.names = FALSE)
  
  #correct format
  df_new <- format(df_new, scientific = scientific)
  
  return(df_new)
}




PlotHist <- function(df_org, df_new, transformation = FALSE, bins = 30){
  #concatenate old and new data
  rbind(cbind(df_org %>% pivot_longer(-1), "Source" = "org"),
        cbind(df_new %>% pivot_longer(-1), "Source" = "new")) %>%
    filter(value > 0) %>%
    #transform data if required
    rowwise() %>%
    mutate(value = ifelse(transformation == FALSE, value, log10(value))) %>%
    #bin data
    group_by(cut(value, bins), Source) %>%
    summarise("number" = n(), .groups = "drop") %>%
    mutate("value" = gsub(".(.*)]", "\\1", .[[1]])) %>%
    separate(value, into = c("from", "to"), sep = ",", convert = TRUE) %>%
    rowwise() %>%
    mutate("mean" = mean(c(from, to), na.rm = TRUE)) %>%
    #correct for differences in data set size
    mutate("number" = ifelse(Source == "new", number*(ncol(df_org)/ncol(df_new)), number)) %>%
    ggplot(aes(mean, number, fill = Source)) +
    geom_col(position = "dodge") +
    theme_minimal()
}
