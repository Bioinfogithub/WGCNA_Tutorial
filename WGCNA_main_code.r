# script to perform WGCNA
setwd("C:\\Users\\vaibhav\\Downloads\\WGCNA of covid dataset")

library(WGCNA) #for weighted gene correlation network analysis 
library(DESeq2) #for differential gene expression analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("devtools")
library(devtools) # tool to make developing  r pacakage easier
ll
#devtools::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot) #CorLevelPlot provides a quick and colourful way to visualise statistically significant correlations between any combination of categorical and continuous variables.
#BiocManager::install("GEOquery") #GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor                                           

library(BiocManager)#to install and manage packages from the Bioconductor project for the statistical analysis and comprehension of high-throughput genomic data.                 
                                                                                                                                                                                                                                  
library(GEOquery)
library(tidyverse) #tidyverse package actually contains other packages (dplyr, ggplot2, etc.) and you’ll see that when you load the tidyverse package using library()                                  
#Tidyverse is a collection of essential R packages for data science

library(CorLevelPlot)
library(gridExtra)#Provides a number of user-level functions to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables

allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data ------------------------------------------------
data <- read.delim('GSE152418_p20047_Study1_RawCounts.txt', header = T)

# get metadata using geo function
geo_id <- "GSE152418"                                                                                                                                           
gse <- getGEO(geo_id, GSEMatrix = TRUE)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

phenoData <- pData(phenoData(gse[[1]]))  #get phenodata                                                                                                                                                                                                                                                               
head(phenoData)                                                      
phenoData <- phenoData[,c(1,2,46:50)]#sub-setting phenodata & get column 1,2, & 46 to 50

# prepare data
data[1:10,1:10] #looking at first 10 rows & first 10 columns

data <- data %>% 
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>% #all sample in col sample and count into col and exclude 1st column (ensemble id)
  mutate(samples = gsub('\\.', '-', samples)) %>% # change sample name
  inner_join(., phenoData, by = c('samples' = 'title')) %>% # merge both the dataframes by col
  select(1,3,4) %>% # keep 1,3,4 column
  spread(key = 'geo_accession', value = 'counts') %>% #column geo_accession
  column_to_rownames(var = 'ENSEMBLID')

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK# checking the outliers(FALSE= we have outliers, TRUE= we don't have outliers)

table(gsg$goodGenes)#checking the outliers & good gene                                                                                                                                                                                                                                                                                                                                                                                                       
table(gsg$goodSamples)#checking the samples (TRUE means good sample)

# remove genes that are detected as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2 #another methods of detecting outliers

pca <- prcomp(t(data))                                                                                                                                                                                                           
pca.dat <- pca$x

pca.var <- pca$sdev^2 #variance
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2) # percentage of variance explained and rounds upto 2 digits

pca.dat <- as.data.frame(pca.dat)# data frame

# plotting the pca data
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %')) # three samples are outliers in plot

###################################### method for Saving the plot in terminal  #################################
#Specify files to save your image using a function such as jpeg(), png(), svg() or pdf(). Additional argument indicating the width and the height of the image can be also used.
#Create the plot
#Close the file with dev.off()

# Open a pdf file
pdf("rplot.pdf") 
# 2. Create a plot
plot(x = my_data$wt, y = my_data$mpg,
     pch = 16, frame = FALSE,
     xlab = "wt", ylab = "mpg", col = "#2E9FDF")
# Close the pdf file
dev.off()


Or use this:

# 1. Open jpeg file
jpeg("rplot.jpg", width = 350, height = "350")
# 2. Create the plot
plot(x = my_data$wt, y = my_data$mpg,
     pch = 16, frame = FALSE,
     xlab = "wt", ylab = "mpg", col = "#2E9FDF")
# 3. Close the file
dev.off()

####################################################################################################################

### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples(3 samples)
samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization ----------------------------------------------------------------------
#if we have micro array data or rnaseq data (in rpkm, fpkm) we need to log transform them

# create a deseq2 data set

# exclude outlier samples in phenodata
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)

  
# fixing column names in col Data
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))# removing additional characters & spaces from names
names(colData) <- gsub('\\s', '_', names(colData))# replacing spaces from col names

# making the row names and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))# checking they are in same order


# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not specifying model



## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNA seq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]# 24 samples
nrow(dds75) # 13284 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t() # transforming the data(rows to be samples, & columns should be gene-id)


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

#plot for mean-connectivity
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
  

#after looking onto plot we have picked out soft -threshold and now we need to perform various steps to identify modules

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor #correlation function
cor <- WGCNA::cor


# memory estimate w.r.t blocksize(blockwiseModules= default, but according to our data we need to set the values)
#default value = 5000,4gb=8-10000,  16gb RAM= 20,000` genes` ,& 232gigs= 30,000 genes
#i set 14000 because i want to  proceesed all the genes in one block instead of that being splitted into multiple blocks

bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,#for reprodeucibility
                 verbose = 3)


cor <- temp_cor #coorelation into new variable


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes) #module-eigen name with color


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),# labels
                    dendroLabels = FALSE,
                    addGuide = TRUE, #plotting parameters
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module





# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% #covid19 infected have value=1,rest are 0
  select(8)


# binarize categorical variables

colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))

severity.out <- binarizeCategoricalColumns(colData$severity,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)


traits <- cbind(traits, severity.out)# combine data with traits data


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')#calculate the coorelation b/w module_eigengenes& traits using pearson coorelation
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples) # p value calculation



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')# combine data from module_eigengenes & traits by rownames

head(heatmap.data)# LOOK AT THE DATA

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')#convert  columns to rownames



#HEATMAP PLOT
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22],# we needs to adjust the number of columns
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))#define colors
#after looking at the plot the level of sigificance is indicated by the no.of astrics(***= means high significance associated with disease, severity, icu)

#which gene belong to what module
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% # to get the gene that belong to module turquoise
  filter(`bwnet$colors` == 'turquoise') %>% #using this genes we can perform further analysis
  rownames()



# 6B. Intramodular analysis: Identifying driver genes ---------------



# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')#module.membership.measures as pearson correlation
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples) # p-values for these measures


module.membership.measure.pvals[1:10,1:10] #memberships and p-value


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')# coorelating expression data with trait of interest, cor method as pearson
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples) # p-values for these gene significance corelations


gene.signf.corr.pvals %>% #top genes that are significantly associated with patients having sever covid19
  as.data.frame() %>% 
  arrange(V1) %>% #(columns is p-values)arrange the p-valus in ascending order
  head(25)# top 25 significant genes

#we can extract these genes ids and convert them into gene symbols & further look into these genes 
# and perform downstream analysis using these genes.

# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.





