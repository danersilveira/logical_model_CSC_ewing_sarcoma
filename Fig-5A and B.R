
library(singscore)
library(msigdbr)
library(msigdb)
library(tidyverse)
library(GSEABase)



# Input Data
# Please check working_data directory for information about how to download the normalized TPM matrix and metadata. See SKNMC-RNA-seq.txt
GSE224004_norm_filt <- readRDS('working_data/GSE224004_tpm_filt.rds')
metadata_GSE224004_filt <- readRDS('working_data/metadata_GSE224004_filt.rds')



###########################################    Signature EWS:FLI1



c2_df <- msigdbr(
  species = "Homo sapiens",
  collection = "C2"
)


ews_kinsey <- c2_df %>%
  filter(gs_name %in% c("KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN", "KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_UP"))

# separar em listas de genes (formato comum para GSEA)
ews_kinsey_list <- ews_kinsey %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)




# Kinsey_UP signature
Kinsey_UP <- ews_kinsey_list$KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_UP[ews_kinsey_list$KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_UP != '']
Kinsey_UP <- unique(Kinsey_UP)
Kinsey_UP

# create GeneSet object
Kinsey_UP_sig <- GeneSet(as.character(Kinsey_UP), setName = 'Kinsey_UP', geneIdType = SymbolIdentifier())
Kinsey_UP_sig

# Kinsey_UP signature
Kinsey_DN <- ews_kinsey_list$KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN[ews_kinsey_list$KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN != '']
Kinsey_DN <- unique(Kinsey_DN)
Kinsey_DN

# create GeneSet object
Kinsey_DN_sig <- GeneSet(as.character(Kinsey_DN), setName = 'Kinsey_DN', geneIdType = SymbolIdentifier())
Kinsey_DN_sig


#rank genes based on expression (logRPKM)
eranks = rankGenes(GSE224004_norm_filt)

#compute epithelial scores
Kinsey_UP_score = simpleScore(eranks, Kinsey_UP_sig)
head(Kinsey_UP_score)

#compute mesenchymal scores
Kinsey_DN_score = simpleScore(eranks, Kinsey_DN_sig)
head(Kinsey_DN_score)

#groups
metadata_GSE224004_filt <- metadata_GSE224004_filt[rownames(metadata_GSE224004_filt) != 'RNAseq_SKNMC_shScramble_48h_rep3', ]
grupo <- paste(metadata_GSE224004_filt$treatment, metadata_GSE224004_filt$time, sep = '-')
grupo


#load MSigDB gene-sets
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM')
hemt_sig = msigdb.hs[['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']]
hemt_score = simpleScore(eranks, hemt_sig)


## GOBP_MESENCHYMAL_STEM_CELL_DIFFERENTIATION 
M_stem_diff_sig = msigdb.hs[['GOBP_MESENCHYMAL_CELL_DIFFERENTIATION']]


M_stem_diff_score = simpleScore(eranks, M_stem_diff_sig)


# Formatting
scores_A <- M_stem_diff_score %>% tibble::rownames_to_column("sample")
scores_B <- Kinsey_UP_score %>% tibble::rownames_to_column("sample")
scores_C <- Kinsey_DN_score %>% tibble::rownames_to_column("sample")

# Merging
scores_combined <- scores_A %>%
  inner_join(scores_B, by = "sample", suffix = c("_A", "_B")) %>%
  inner_join(scores_C, by = "sample", suffix = c("_A", "_C"))

scores_combined$grupo <- grupo









