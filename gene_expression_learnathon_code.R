source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)

gds <- getGEO("GDS5056")  # 3 obese and 3 non-obese patient's adipose stem cells
str(gds)

gds <- getGEO("GDS5167") # 6 obese diabetic and 6 obese non-diabetic visceral adipose tissue CD14+ cells
str(gds)

gds <- getGEO("GDS6177") # 5 patients- acute alcohol consumption effect on whole blood (time series) control group
str(gds)
gds <- getGEO("GDS4938") # 5 patients- acute alcohol consumption effect on whole blood (time series) alcohol group
str(gds)

gds <- getGEO("GDS4758") # alzheimer's decesead brains 78 patients
str(gds)

gds <- getGEO("GDS5363") # osteoarthritis onset- 139 samples
str(gds)

gds <- getGEO("GDS5037") # asthma development 108 samples (control, mild, severe)
str(gds)

gds <- getGEO("GDS4971") # sepsis survivors and non-survivors- time series 163 samples
str(gds)

gds <- getGEO("GDS4274") # 130 samples, pediatric septic shock clinical subgroups
str(gds)

gds <- getGEO("GDS5393") # bipolar disorder and response to lithium treatment, 120 samples
str(gds)


#EXAMPLE DATASET FROM HERE





install.packages("rpart")
library(rpart)
myTree=rpart(TrainingLabels~.,cbind(TrainingData[,1:10000],TrainingLabels),method="class")
plot(myTree)
text(myTree)

#Installing necessary packages to run rattle
install.packages("RGtk2")
install.packages("rattle")
install.packages("rpart.plot")

#Load in the packages to create a fancified version of your tree
library(rattle)
library(rpart.plot)
library(RColorBrewer)

fancyRpartPlot(myTree)


install.packages("randomForest")
library(randomForest)
fit=randomForest(TrainingData,y=TrainingLabels,ntree=500)
prediction=predict(fit, TestingData,type = "response")
table(prediction,TestingLabels)


#k means example
data(iris)
summary(iris)
k<-kmeans(iris[,1:4],centers=3)
plot(iris[,1], iris[,2], col=as.factor(k$cluster)) 
plot(iris[,1], iris[,2], col=as.factor(iris[,5]))

plot(iris[,1:4], col=as.factor(iris[,5]))
plot(iris[,1:4], col=as.factor(k$cluster))


# k nearest neighbors
install.packages('class')
library(class)
prediction=knn(TrainingData, TestingData, TrainingLabels, k = 5)
table(prediction,TestingLabels)





# differential gene expression p-values

pvals=apply(temp,MARGIN = 2,function(x){t.test(x[1:64],x[65:122])$p.value})
pvals_wilcox=apply(temp,MARGIN = 2,function(x){wilcox.test(x[1:64],x[65:122])$p.value})
geneNames=table[names(pvals),2]

adjusted_pvals=pvals*length(geneNames)
adjusted_pvals_wilcox=pvals_wilcox*length(geneNames)

fdr_vals=pvals*length(geneNames)/length(which(pvals<0.05))
fdr_vals_wilcox=pvals_wilcox*length(geneNames)/length(which(pvals_wilcox<0.05))

topGenes=geneNames[which(adjusted_pvals<0.05)]
topGenes=geneNames[which(adjusted_pvals_wilcox<0.05)]

topGenes=geneNames[which(fdr_vals<0.05)]
topGenes=geneNames[which(fdr_vals_wilcox<0.05)]




#gene ontology terms

source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
biocLite("org.Hs.eg.db")
library("GO.db")
library("org.Hs.eg.db")
#The package GO.db contains organism-independent information about Gene Ontology terms from geneontology.org, while org.Hs.eg.db contains specific information about human genes.

keytypes(GO.db)
select(GO.db, keytype="TERM", keys="galactose metabolic process", columns = "GOID")
select(GO.db, keys="GO:0006012", columns=c("TERM","DEFINITION"), keytype="GOID")
keytypes(org.Hs.eg.db)
keys(org.Hs.eg.db,keytype="ALIAS")
columns(org.Hs.eg.db)

GN.table = select(org.Hs.eg.db, keys="GO:0006012", columns="ALIAS", keytype="GO")
all.annotated = keys(org.Hs.eg.db,keytype="ALIAS")
