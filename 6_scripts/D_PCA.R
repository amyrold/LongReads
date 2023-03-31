#PCA stands for Principal component Analysis
#PCA plots are a dimensional reduction method in which the dimensions of a data set are reduced. 
#You essentially reduce the number of variables of a data set while preserving as much information as possible. 
#You’re trading some accuracy for more simplicity and a better visual representation of the correlation. 

install.packages("ggfortify") #package used to create PCA plots
library(ggfortify)

df <- iris [1:4] #storing the iris data set into a data frame, a tabular data structure
#The iris data set is a data set of multiple measurements and taxonomic names of 3 different plant species
df


pca_res <- prcomp(df, scale. = TRUE) #this scales the iris data by dividing its values by std. dev and adds PCA labels
pca_res 
autoplot(pca_res) #Creates a basic PCA plot

#PCA plots can be customized and edited in many ways shown below
autoplot(pca_res, data = iris, colour = "Species") #colors points based on species
autoplot(pca_res, data = iris, colour = "Species", label = TRUE, label.size = 3) #draws labels
autoplot(pca_res, data = iris, colour = "Species", shape = FALSE, label.size = 3) #removes points
autoplot(pca_res, data = iris, colour = "Species", loadings = TRUE) #draws eigenvectors
autoplot(pca_res, data = iris, colour = "Species", loadings = TRUE, loadings.color = 'blue',
         loadings.label = TRUE, loadings.size = 3) #labels eigenvectors


library(cluster) #clustering package
autoplot(clara(iris[-5], 3)) #allows clustering of data points
autoplot(fanny(iris[-5], 3), frame = TRUE) #Outlines each cluster
autoplot(pam(iris[-5], 3), frame = TRUE, frame.type = 'norm') #Another way to outline each cluster


#A Pseudo Data set of Edit distances of Intragenomic fasta files and Intergenomic fasta files
#Below is code that plots the edit distance data into a PCA graph and clusters them
Intragenomic = c(14, 20, 6, 9 , 10, 8, 21, 5, 16, 17, 35, 37, 43, 59, 58, 29, 28, 27, 42, 61, 78, 70)
Intergenomic = c(23, 24, 1, 11, 12, 19, 18, 26, 13, 22, 66, 42, 41, 59, 74, 80, 81, 50, 53, 63, 40, 57)
fasta_frame <- data.frame(Intragenomic, Intergenomic) #storing edit distance values as a data frame
fasta_frame[1:2]

pca_res <- prcomp(fasta_frame, scale. = TRUE) #scaling edit distance data
autoplot(clara(fasta_frame[-4], 2)) #basic PCA plot of edit distance data
autoplot(fanny(fasta_frame[-4], 2), frame = TRUE) ##basic clustering
