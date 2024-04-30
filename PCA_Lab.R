
#Reading the data 
data <- read.csv("C:/Users/DELL/Downloads/data.csv", row.names = 1)

str(data)
View(data)

#load the library
library("factoextra")
library("FactoMineR")

#standardize the values then creates a new principal component table
pca.data <- PCA(data[,-1], scale.unit = TRUE, graph = FALSE)

#we need to use the fviz_eig() function.
fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))

# understand the correlation between the samples
fviz_pca_var(pca.data, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

#start by plotting the cell types 
pca.data <- PCA(t(data[,-1]), scale.unit = TRUE, graph = FALSE)

#visualization
fviz_pca_ind(pca.data, col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
             repel = TRUE)

#load the library
library(ggpubr) 


a <- fviz_pca_ind(pca.data, col.ind = "cos2", 
                  gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                  repel = TRUE)

#we can use ggpar() 
ggpar(a,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())

#plot the gene
pca.data <- PCA(data[,-1], scale.unit = TRUE,ncp = 2, graph = FALSE)

#color the gene
data$lineage <- as.factor(data$lineage)


#load the library
library(RColorBrewer)

nb.cols <- 3
mycolors <- colorRampPalette(brewer.pal(3, "Set1"))(nb.cols)

#Lineage column
a <- fviz_pca_ind(pca.data, col.ind = data$lineage,
                  palette = mycolors, addEllipses = TRUE)

#we will use the ggpar() function to add labels
ggpar(a,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cell type", legend.position = "top",
      ggtheme = theme_minimal())





#Q1. After you run PCA() and fviz_eig() functions, you will see the plot. What is the number of the first two principal component variations? 
  
# Answer : The first two principal components explain 81.7% of the variation.

#Q2. After you follow the tutorial, you will get the gene PCA plot. Based on the gene expression plot discussed above, identify the marker genes that are highly expressed in trophoectoderm (TE) cells, epiblast (EPI) cells, and primitive endoderm (PE) cells. Explain why these genes are significant in identifying specific cell types. 

#Answer : we can see that there are some genes are associated with a specific type of cell. These genes could act as marker to identify these cell types. In TE cells the KRT18, KRT8 and S100A16 are highly expresses compare to the other genes. In EPI cells the DPPA5, IFITM1 , MT1X and UPP1 are highly expressed compare to other genes. And in PE cells only APOA1 is highly expressed compare to the other genes.
