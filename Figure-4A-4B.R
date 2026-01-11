# Read packages
library(dplyr)
library(limma)
library(fgsea)
library(msigdbr)


# Input Data
# Please check working_data directory for information about how to download and process the rds files below. Two txt files called Ewing_cell_lines_microarray.txt and Metadata_GSE176190.txt 
# are provided related to Ewing_cell_lines_microarray.rds and Metadata_GSE176190-90-cells.rds, respectively.
GSE176190 <- readRDS('working_data/Ewing_cell_lines_microarray.rds') 
metadata_GSE176190 <- readRDS('Metadata_GSE176190-90-cells.rds')


################ limma 

# Test correspondence between GSE176190 colnames and metadata
all(colnames(GSE176190) == metadata_GSE176190$geo_accession)

# If TRUE
colnames(GSE176190) <- metadata_GSE176190$title


# Defining groups
Groups <- unlist(metadata_GSE176190$condition)

# design matrix
design <- model.matrix(~factor(Groups) -1)

# Defining colnames
colnames(design) <- levels(as.factor(Groups))

# Fit
fit <- lmFit(GSE176190, design)

colnames(fit)

# Making contrast
constrast.matrix <- makeContrasts('A673_EETS_low-A673_EETS_high', 'CHLA10_EETS_low-CHLA10_EETS_high','EW1_EETS_low-EW1_EETS_high', 'EW22_EETS_low-EW22_EETS_high',
                                  'EW24_EETS_low-EW24_EETS_high', 'EW7_EETS_low-EW7_EETS_high', 'MHHES1_EETS_low-MHHES1_EETS_high', 'MIC_EETS_low-MIC_EETS_high',
                                  'POE_EETS_low-POE_EETS_high', 'RDES_EETS_low-RDES_EETS_high', 'RH1_EETS_low-RH1_EETS_high','SKES1_EETS_low-SKES1_EETS_high',
                                  'SKNMC_EETS_low-SKNMC_EETS_high', 'TC32_EETS_low-TC32_EETS_high', 'TC71_EETS_low-TC71_EETS_high',
                                  levels = design)
fitc <- contrasts.fit(fit, constrast.matrix)
fitc <- eBayes(fitc)
res <- topTable(fitc, number=Inf)

summary(decideTests(fitc, method="global"))

# Filtering by cell line
res_A673 <-  topTable(fitc, number=Inf, coef=1)
res_CHLA10 <-  topTable(fitc, number=Inf, coef=2)
res_EW1 <-  topTable(fitc, number=Inf, coef=3)
res_EW22 <- topTable(fitc, number=Inf, coef=4)
res_EW24 <-  topTable(fitc, number=Inf, coef=5)
res_EW7 <-  topTable(fitc, number=Inf, coef=6)
res_MHHES1 <-  topTable(fitc, number=Inf, coef=7)
res_MIC <- topTable(fitc, number=Inf, coef=8)
res_POE <-  topTable(fitc, number=Inf, coef=9)
res_RDES <-  topTable(fitc, number=Inf, coef=10)
res_RH1 <-  topTable(fitc, number=Inf, coef=11)
res_SKES1 <-  topTable(fitc, number=Inf, coef=12)
res_SKNMC <-  topTable(fitc, number=Inf, coef=13)
res_TC32 <-  topTable(fitc, number=Inf, coef=14)
res_TC71 <-  topTable(fitc, number=Inf, coef=15)

###################################################  Enrichment  ###############################################################################



# Set seed
set.seed(21)

# Download Hallmark
hallmark_df <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

# Converting to list
hallmakr_list <- hallmark_df %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Download Gene ontology pathways
ontology_df <- msigdbr(
  species  = "Homo sapiens",
  collection = "C5"
)


pathways <- c('GOBP_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION', 'GOBP_REGULATION_OF_STEM_CELL_DIFFERENTIATION','GOBP_STEM_CELL_DIFFERENTIATION',
              'GOBP_REGULATION_OF_MESENCHYMAL_STEM_CELL_DIFFERENTIATION', 'GOBP_EPITHELIAL_MESENCHYMAL_CELL_SIGNALING', 'GOBP_NEGATIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION', 
              'GOBP_POSITIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION', 'GOBP_MESENCHYMAL_STEM_CELL_DIFFERENTIATION', 'GOBP_MESENCHYMAL_CELL_DEVELOPMENT', 'GOBP_EPITHELIAL_CELL_CELL_ADHESION', 'GOBP_EPITHELIAL_CELL_DEVELOPMENT', 'GOBP_EPITHELIAL_CELL_DIFFERENTIATION',
              'GOBP_MESENCHYMAL_CELL_DIFFERENTIATION', 'GOBP_MESENCHYMAL_CELL_DEVELOPMENT', 'GOBP_EPITHELIAL_CELL_CELL_ADHESION', 'GOBP_EPITHELIAL_CELL_DEVELOPMENT', 'GOBP_EPITHELIAL_CELL_DIFFERENTIATION',
              'GOBP_MESENCHYMAL_CELL_DIFFERENTIATION')


# Filtering
ontology_df_fit <- ontology_df[ontology_df$gs_name %in% pathways, ]

# Converting to list
ontology_list <- ontology_df_fit %>%
  split(x = .$gene_symbol, f = .$gs_name)



########################################################################## A673 #######################################################



stats_A673 <- res_A673$logFC
names(stats_A673) <- rownames(res_A673)
stats_A673 <- stats_A673[is.finite(stats_A673)]
stats_A673 <- sort(stats_A673, decreasing = TRUE)

# hallmark

fgsea_res_hall_A673 <- fgsea(pathways = hallmakr_list, 
                   stats = stats_A673, 
                   eps = 0.0,
                   minSize = 15,
                   maxSize = 500)



# GO BP STEM


fgsea_res_GO_A673 <- fgsea(pathways = ontology_list, 
                             stats = stats_A673, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)



########################################################################## CHLA10 #######################################################



stats_CHLA10 <- res_CHLA10$logFC
names(stats_CHLA10) <- rownames(res_CHLA10)
stats_CHLA10 <- stats_CHLA10[is.finite(stats_CHLA10)]
stats_CHLA10 <- sort(stats_CHLA10, decreasing = TRUE)

# hallmark

fgsea_res_hall_CHLA10 <- fgsea(pathways = hallmakr_list, 
                             stats = stats_CHLA10, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)


# GO BP STEM


fgsea_res_GO_CHLA10 <- fgsea(pathways = ontology_list, 
                           stats = stats_CHLA10, 
                           eps = 0.0,
                           minSize = 15,
                           maxSize = 500)


########################################################################## EW1 #######################################################



stats_EW1 <- res_EW1$logFC
names(stats_EW1) <- rownames(res_EW1)
stats_EW1 <- stats_EW1[is.finite(stats_EW1)]
stats_EW1 <- sort(stats_EW1, decreasing = TRUE)

# hallmark

fgsea_res_hall_EW1 <- fgsea(pathways = hallmakr_list, 
                               stats = stats_EW1, 
                               eps = 0.0,
                               minSize = 15,
                               maxSize = 500)



# GO BP STEM


fgsea_res_GO_EW1 <- fgsea(pathways = ontology_list, 
                             stats = stats_EW1, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)




########################################################################## EW22 #######################################################



stats_EW22 <- res_EW22$logFC
names(stats_EW22) <- rownames(res_EW22)
stats_EW22 <- stats_EW22[is.finite(stats_EW22)]
stats_EW22 <- sort(stats_EW22, decreasing = TRUE)

# hallmark

fgsea_res_hall_EW22 <- fgsea(pathways = hallmakr_list, 
                            stats = stats_EW22, 
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)


# GO BP STEM


fgsea_res_GO_EW22 <- fgsea(pathways = ontology_list, 
                          stats = stats_EW22, 
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)




########################################################################## EW24 #######################################################



stats_EW24 <- res_EW24$logFC
names(stats_EW24) <- rownames(res_EW24)
stats_EW24 <- stats_EW24[is.finite(stats_EW24)]
stats_EW24 <- sort(stats_EW24, decreasing = TRUE)

# hallmark

fgsea_res_hall_EW24 <- fgsea(pathways = hallmakr_list, 
                             stats = stats_EW24, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)



# GO BP STEM


fgsea_res_GO_EW24 <- fgsea(pathways = ontology_list, 
                           stats = stats_EW24, 
                           eps = 0.0,
                           minSize = 15,
                           maxSize = 500)






########################################################################## EW7 #######################################################



stats_EW7 <- res_EW7$logFC
names(stats_EW7) <- rownames(res_EW7)
stats_EW7 <- stats_EW7[is.finite(stats_EW7)]
stats_EW7 <- sort(stats_EW7, decreasing = TRUE)

# hallmark

fgsea_res_hall_EW7 <- fgsea(pathways = hallmakr_list, 
                             stats = stats_EW7, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)

df_result <- fgsea_res_hall_EW7
df_result <- df_result[order(df_result$padj), ]

top_pathways_hall_EW7 <- df_result %>% 
  filter(padj < 0.05) %>%
  arrange(padj) %>% 
  head(10) 

# GO BP STEM


fgsea_res_GO_EW7 <- fgsea(pathways = ontology_list, 
                           stats = stats_EW7, 
                           eps = 0.0,
                           minSize = 15,
                           maxSize = 500)




########################################################################## MHHES1 #######################################################



stats_MHHES1 <- res_MHHES1$logFC
names(stats_MHHES1) <- rownames(res_MHHES1)
stats_MHHES1 <- stats_EW7[is.finite(stats_MHHES1)]
stats_MHHES1 <- sort(stats_MHHES1, decreasing = TRUE)

# hallmark

fgsea_res_hall_MHHES1 <- fgsea(pathways = hallmakr_list, 
                            stats = stats_MHHES1, 
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)



# GO BP STEM


fgsea_res_GO_MHHES1 <- fgsea(pathways = ontology_list, 
                          stats = stats_MHHES1, 
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)




########################################################################## MIC #######################################################



stats_MIC <- res_MIC$logFC
names(stats_MIC) <- rownames(res_MIC)
stats_MIC <- stats_EW7[is.finite(stats_MIC)]
stats_MIC <- sort(stats_MIC, decreasing = TRUE)

# hallmark

fgsea_res_hall_MIC <- fgsea(pathways = hallmakr_list, 
                               stats = stats_MIC, 
                               eps = 0.0,
                               minSize = 15,
                               maxSize = 500)



# GO BP STEM


fgsea_res_GO_MIC <- fgsea(pathways = ontology_list, 
                             stats = stats_MIC, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)


########################################################################## POE #######################################################



stats_POE <- res_POE$logFC
names(stats_POE) <- rownames(res_POE)
stats_POE <- stats_EW7[is.finite(stats_POE)]
stats_POE <- sort(stats_POE, decreasing = TRUE)

# hallmark

fgsea_res_hall_POE <- fgsea(pathways = hallmakr_list, 
                            stats = stats_POE, 
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)



# GO BP STEM


fgsea_res_GO_POE <- fgsea(pathways = ontology_list, 
                          stats = stats_POE, 
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)


########################################################################## RDES #######################################################



stats_RDES <- res_RDES$logFC
names(stats_RDES) <- rownames(res_RDES)
stats_RDES <- stats_EW7[is.finite(stats_RDES)]
stats_RDES<- sort(stats_RDES, decreasing = TRUE)

# hallmark

fgsea_res_hall_RDES <- fgsea(pathways = hallmakr_list, 
                            stats = stats_RDES, 
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)



# GO BP STEM


fgsea_res_GO_RDES <- fgsea(pathways = ontology_list, 
                          stats = stats_RDES, 
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)




########################################################################## RH1 #######################################################



stats_RH1 <- res_RH1$logFC
names(stats_RH1) <- rownames(res_RH1)
stats_RH1 <- stats_EW7[is.finite(stats_RH1)]
stats_RH1 <- sort(stats_RH1, decreasing = TRUE)

# hallmark

fgsea_res_hall_RH1 <- fgsea(pathways = hallmakr_list, 
                             stats = stats_RH1, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)



# GO BP STEM


fgsea_res_GO_RH1 <- fgsea(pathways = ontology_list, 
                           stats = stats_RH1, 
                           eps = 0.0,
                           minSize = 15,
                           maxSize = 500)


########################################################################## SKES1 #######################################################



stats_SKES1 <- res_SKES1$logFC
names(stats_SKES1) <- rownames(res_SKES1)
stats_SKES1 <- stats_SKES1[is.finite(stats_SKES1)]
stats_SKES1 <- sort(stats_SKES1, decreasing = TRUE)

# hallmark

fgsea_res_hall_SKES1 <- fgsea(pathways = hallmakr_list, 
                            stats = stats_SKES1, 
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)



# GO BP STEM


fgsea_res_GO_SKES1 <- fgsea(pathways = ontology_list, 
                          stats = stats_SKES1, 
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)



########################################################################## SKNMC #######################################################



stats_SKNMC <- res_SKNMC$logFC
names(stats_SKNMC) <- rownames(res_SKNMC)
stats_SKNMC <- stats_SKNMC[is.finite(stats_SKNMC)]
stats_SKNMC <- sort(stats_SKNMC, decreasing = TRUE)

# hallmark

fgsea_res_hall_SKNMC <- fgsea(pathways = hallmakr_list, 
                              stats = stats_SKNMC, 
                              eps = 0.0,
                              minSize = 15,
                              maxSize = 500)



# GO BP STEM


fgsea_res_GO_SKNMC <- fgsea(pathways = ontology_list, 
                            stats = stats_SKNMC, 
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)


########################################################################## TC32 #######################################################



stats_TC32<- res_TC32$logFC
names(stats_TC32) <- rownames(res_TC32)
stats_TC32 <- stats_TC32[is.finite(stats_TC32)]
stats_TC32 <- sort(stats_TC32, decreasing = TRUE)

# hallmark

fgsea_res_hall_TC32 <- fgsea(pathways = hallmakr_list, 
                              stats = stats_TC32, 
                              eps = 0.0,
                              minSize = 15,
                              maxSize = 500)



# GO BP STEM


fgsea_res_GO_TC32 <- fgsea(pathways = ontology_list, 
                            stats = stats_TC32, 
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)





########################################################################## TC71 #######################################################



stats_TC71 <- res_TC71$logFC
names(stats_TC71) <- rownames(res_TC71)
stats_TC71 <- stats_TC71[is.finite(stats_TC71)]
stats_TC71 <- sort(stats_TC71, decreasing = TRUE)

# hallmark

fgsea_res_hall_TC71 <- fgsea(pathways = hallmakr_list, 
                             stats = stats_TC71, 
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)



# GO BP STEM


fgsea_res_GO_TC71 <- fgsea(pathways = ontology_list, 
                           stats = stats_TC71, 
                           eps = 0.0,
                           minSize = 15,
                           maxSize = 500)

