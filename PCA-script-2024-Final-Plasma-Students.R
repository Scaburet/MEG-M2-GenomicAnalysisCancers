

##############################
# 1. Packages and workspace  #
##############################

# Let's import the packages we need !

# 1. Which packages are installed on PLASMA server?
# Is FactoMineR present into the library?
# NB : if needed, to close this list of librairies, and go back to the command prompt: q


## Which packages are already available in your session?


## Import FactoMineR package into your workspace


## Check that FactoMineR has been imported into your R workspace



####################################
# 2. Choose your working directory #
####################################

## Write the command to know what your actual working directory is:


# on adenine, the absolute path should be like this "/srv/home/login/meg_m2_gac" 
# (login being the one that was provided).
# alternatively, you can use the symbol "~" which stands for your home:
#setwd("~/meg_m2_gac")


###########################################################################
# Create a dataframe in your working space that contains data from a file #
###########################################################################

### To import data: 
##		- from a file in a table format
##		- to a data frame


## read.table function
# Visit the following URL to understand the following command line:
# https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/read.table

# Then, execute the following command:
accdata <- read.table("6-GSE10927-norm.txt", header=TRUE, sep="\t", dec=".", row.names=1)



#########################################################################
# Perform the Principal Component Analysis using the FactoMineR package #
#########################################################################

## PCA function has been imported into your workspace with the FactoMineR librairy
# Go to the following URL to understand the argumant provided in the following function:
# https://cran.r-project.org/web/packages/FactoMineR/FactoMineR.pdf

res.pca <- PCA(accdata[,1:1905], quali.sup = 1901:1904, quanti.sup = 1905, graph = FALSE)



#############################################################
################## Create plots of the PCA ##################
#############################################################


# Plot of individuals
######################


## plot function (and arguments between brackets)
# Let's plot the PCA result for individuals ("ind" argument) !

plot(res.pca, choix = "ind", cex=0.6, invisible = "quali")

# invisible = "quali" indicates that qualitative variables centroid shall not be drawn on the plot

## To save your plot:
# Informations provided at https://cran.r-project.org/web/packages/jpeg/jpeg.pdf

jpeg("PCA_ind.jpg", width = 700, height = 700)
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali")
dev.off()


# Plot of individuals with annotations (supplementary variables)
# Let's try to find a biological meaning to this PCA !!!!!
#################################################################################################################

# Add the information of supplementary qualitative variables
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", autoLab = "y")
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "clusterforACCs", autoLab = "y")
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Side", autoLab = "y")

# res.pca is the name of the object used by the plot function
# invisible = "quali" indicates that qualitative variables centroid shall not be drawn on the plot
# habillage colors individuals among a categorical variable (here "Diagnosis", which is the name of the column)

# Add the information of supplementary qualitative variables, for axes 2 and 3
plot(res.pca, axes=c(2,3), choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", autoLab = "y")
plot(res.pca, axes=c(2,3), choix = "ind", cex=0.6, invisible = "quali", habillage = "clusterforACCs", autoLab = "y")
plot(res.pca, axes=c(2,3), choix = "ind", cex=0.6, invisible = "quali", habillage = "Side", autoLab = "y")

# res.pca is the name of the object used by the plot function
# invisible = "quali" indicates that qualitative variables centroid shall not be drawn on the plot
# habillage colors individuals among a categorical variable (here "Diagnosis", which is the name of the column)
# axes=c(2,3) plot the 2nd and 3rd maximum inertia axes


# Try to find qualitative variables represented by axes:
##########################################################

# To quantify the correlation between qualitative variables and coordinates on the axes:
# NB: Cos2 value representes the quality of the representation of the qualitative variable by the axis
# close to 1: good representation, close to zero: bad representation
quali_results <- res.pca$quali
print(quali_results$cos2)


# Plot of variables, try to find which mRNA levels are correlated to biological information
###########################################################################################

plot(res.pca, axes = c(1,2), cex=0.6, choix="var", select="contrib 30", unselect = 1)

# select="contrib 30" select the first 30 variables that contribute to PC1 and PC2 (axes illustrated by default)

# Have a look on mRNA contributing to axes 2 and 3
plot(res.pca, axes = c(2,3), choix="var", cex=0.6, select="contrib 15", unselect = 1)



# To print the variables (mRNA) that are well represented by PC1 and PC2:
cos2 <- res.pca$var$cos2
selected_variables_PC1_PC2 <- rownames(cos2)[cos2[, 1] > 0.8 | cos2[, 2] > 0.8]  # Pour PC1 ou PC2
print(selected_variables_PC1_PC2)


#################################################################################
########### Export to files informations important for interpretation ###########
#################################################################################

# To know the informations contained in the res.pca object
res.pca

# Principle Components (coordinates of individual on axes)
write.infile(res.pca$ind$coord, file="PrincipalComponentsInd.txt", sep="\t")

# Correlation of mRNA levels with axes
write.infile(res.pca$var$cor, file="CorrelationVar.txt", sep="\t")

# Coordinates of mRNA levels with
write.infile(res.pca$var$coord, file="CoordinVar.txt", sep="\t")



##############################################################################################
########### Quality of the interpertation depends on your knowledge in the field ! ###########
##############################################################################################


#######################################################################
########### Contribution of Principal Components to Inertia ###########
#######################################################################

## Eigenvalues are the variances associated to each principal components.
## The following eigen values matrix contains the inertia value explained by each of the 64 PCs 
eig.val <- res.pca$eig

## Let us draw a plot
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        ylim = c(0,25),
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")

## Zoom on the first 25 components

jpeg("QualityProj_25C.jpg", width = 700, height = 700)
eig.val <- res.pca$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        ylim = c(0,25),
        xlim = c(0,25),
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")
dev.off()


## Export to a file the eigenvalues
write.infile(res.pca$eig, file="QualityProj_PC.txt", sep="\t")

summary(res.pca, ncp=3, nbelements=Inf, file="pca-essai.txt" )

#######################################
# Identify the 10 samples that contributed the most to the selected axes
# Identify the samples that are correctly represented by the selected axes

plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", select = "contrib 10")
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", select = "cos2 0.2")

# Below, you'll find the same result, without the names of samples
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", select = "contrib 10", label = "none")
plot(res.pca, choix = "ind", cex=0.6, invisible = "quali", habillage = "Diagnosis", select = "cos2 0.2",label = "none")





