library(dplyr)
library(rstatix)


# Input Data
# Please check working_data directory for information about how to download and process the rds files below. One txt file called SKNMC-RNA-seq.txt  
# is provided related to Counts_GSE224004_filt-quantile.rds and metadata_GSE224004_filt-quantile.rds, respectively.
GSE224004_norm_filt <- read.csv('working_data/Counts_GSE224004_filt-quantile.rds',  row.names = 1, header = TRUE)
metadata_GSE224004_filt <- read.csv('working_data/metadata_GSE224004_filt-quantile.rds', sep = ',', row.names = 2)



genes_de_interesse <- c('LIN28B', 'MYC', 'STAT3', 'SOX2', 'POU5F1', 'NT5E', 'PROM1')
Expression_values <- as.matrix(GSE224004_norm_filt[genes_de_interesse, ])


df <- data.frame(
  Gene = rownames(Expression_values),
  Sample = rep(colnames(GSE224004_norm_filt), each = length(genes_de_interesse)),
  Expression = as.numeric(Expression_values),
  Condition = rep(metadata_GSE224004_filt$condition, each = length(genes_de_interesse))
)


# 24h
df_24 <- df %>%
  filter(Condition %in% c("shScramble-24h", "shFLI1-24h")) %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c("shScramble-24h", "shFLI1-24h")
    )
  )



# T-test 
ttest_24 <- df_24 %>%
  group_by(Gene) %>%
  t_test(Expression ~ Condition) %>%   
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") 




# 24h
df_48 <- df %>%
  filter(Condition %in% c("shScramble-48h", "shFLI1-48h")) %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c("shScramble-48h", "shFLI1-48h")
    )
  )



# T-test 
ttest_48 <- df_48 %>%
  group_by(Gene) %>%
  t_test(Expression ~ Condition) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")



