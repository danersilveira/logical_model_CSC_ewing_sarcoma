# Read packages

library(dplyr)
library(rstatix)



# Input Data
# Please check working_data directory for information about how to download and process the rds files below. One txt file called Ewing_cell_lines_microarray.txt 
# is provided related to Ewing_cell_lines_microarray.rds and Metadata_GSE176190-90-cells.rds, respectively.
GSE176190 <- readRDS('working_data/Ewing_cell_lines_microarray.rds') 
metadata_GSE176190 <- readRDS('Metadata_GSE176190-90-cells.rds')



############################################### DEGs analysis ###############################################################################

# Selecting genes of interest
genes_of_interest <- c('LIN28A', 'LIN28B', 'MYC', 'STAT3', 'SOX2', 'POU5F1') 

Expression_values <- as.matrix(GSE176190[genes_of_interest, ])  
df <- data.frame(Gene = rownames(Expression_values),  
                 Sample = rep(colnames(GSE176190), each = length(genes_of_interest)),  
                 Expression = as.numeric(Expression_values), 
                 Condition = rep(metadata_GSE176190$condition, each = length(genes_of_interest)),
                 Cell_line = rep(metadata_GSE176190$`cell.line:ch1`, each = length(genes_of_interest)),
                 Genotype = rep(metadata_GSE176190$`genotype:ch1`, each = length(genes_of_interest)))



# T-test
ttest <- df %>%
  group_by(Cell_line, Gene) %>%
  t_test(Expression ~ Genotype, var.equal = TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")




