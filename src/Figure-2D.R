#Read packages
library("ggplot2")
library("FactoMineR")
library("factoextra")

rm(list=ls(all=TRUE))


#Read data from a file (.csv, .txt, etc...).
modeldata <- read.table('working_data/States_MCA.csv', header = TRUE, sep = ',')

#Check if data were correctly read
head(modeldata)
summary(modeldata)

#MCA package -9 is used to exclude the column PHENOTYPE of the analysis
res.mca <- MCA(modeldata[, -9],  graph=FALSE)

#Coordinates of individuals
ind <- get_mca_ind(res.mca)
ind$coord

#Coordinates of variables
var <- get_mca_var(res.mca)
var$coord

#grp corresponds to the column PHENOTYPE of .csv file
grp <- as.factor(modeldata[, "PHENOTYPE"])

fviz_screeplot(res.mca, addlabels = TRUE)


#Plot 2D graph
plot <- fviz_mca_ind(res.mca,
             #axes = c(1, 1),
             label = "none", # hide individual labels
             labelsize = 9, pointsize = 7,
             habillage = grp, # color by groups
             palette = c("#55308d", "#ffb66c", "#ea7500"),  # Number of colors = Number of phenotypes
             repel = FALSE,
             addEllipses = FALSE, # Concentration ellipses
             ggtheme = theme_minimal(),
             mean.point = FALSE,
) + theme(axis.text = element_text(size = 20)) +  theme(text = element_text(size = 20)) 

plot
