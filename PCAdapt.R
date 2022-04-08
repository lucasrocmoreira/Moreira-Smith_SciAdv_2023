#############################################
#                  PCAdapt                  #
#            Modified on 10-27-20           #
#############################################

#install.packages("pcadapt")
library(pcadapt)

path_to_file <- "SNP-only_simplified.lfmm"
filename <- read.pcadapt(path_to_file, type = "lfmm")

# To choose K, principal component analysis should first be performed with a large enough number of principal components

x <- pcadapt(input = filename, K = 7)

# The 'scree plot' displays in decreasing order the percentage of variance explained by each PC. Up to a constant, it corresponds to the eigenvalues in decreasing order. 
# The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. We recommend to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell's rule).

plot(x, option = "screeplot")

# Let's vizualize the structure

# With integers
poplist.int <- c(rep(1, 10), rep(2, 10), rep(3, 10),rep(4, 10),rep(5, 10),rep(6, 10),rep(7, 10))
# With names
poplist.names <- c(rep("AK",10), rep("MW", 10), rep("NE", 10),rep("NR", 10),rep("NW", 10),rep("SE", 10),rep("SR", 10))
print(poplist.int)
print(poplist.names)

plot(x, option = "scores", pop = poplist.names)

plot(x, option = "scores", i = 1, j =3, pop = poplist.names)

# Looking at population structure beyond K = 4 confirms the results of the scree plot. The fourth and the fifth principal components do not ascertain population structure anymore.

#### COMPUTING THE TEST STATISTIC ###

# For a given SNP, the test statistic is based on the z-scores obtained when regressing SNPs with the K principal components.
# The test statistic for detecting outlier SNPs is the Mahalanobis distance, which is a multi-dimensional approach that measures how distant is a point from the mean.
# Once divided by a constant ?? called the genomic inflation factor, the scaled squared distances D2j/?? should have a chi-square distribution with K degrees of freedom under the assumption that there are no outlier.

# In addition to the number K of principal components to work with, the user can also set the parameter min.maf that corresponds to a threshold of minor allele frequency. By default, the parameter min.maf is set to 5%.

y <- pcadapt(filename, K = 4)

saveRDS(y$pvalues,"p-values.RDS")

# The object x returned by the function pcadapt contains numerical quantities obtained after performing a PCA on the genotype matrix.

################ Visualizing the data ################

hist(y$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# An histogram of p-values confirms that most of the p-values follow an uniform distribution. The excess of small p-values indicates the presence of outliers.

################ Choosing a cutoff for outlier detection ################

# To provide a list of outliers and choose a cutoff for outlier detection, there are several methods that are listed below from the less conservative one to the more conservative one

# For a given ?? (real valued number between 0 and 1), SNPs with q-values less than ?? will be considered as outliers with an expected false discovery rate bounded by ??. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs.

#BiocManager::install("qvalue")
library(qvalue)

qval <- qvalue(y$pvalues)$qvalues
alpha <- 0.01 # expected false discovery rate lower than 1%
outliers <- which(qval < alpha)
length(outliers)
saveRDS(outliers,"Outliers.RDS")

# Create color vector to highlight outliers
length(y$pvalues)
isoutlier <- rep("black",7009778)
isoutlier[outliers] <- "red"

library(dplyr)
pvalues <- data.frame(p=y$pvalues,color=isoutlier)
SNP_reduced <- pvalues %>%
  filter(-log10(pvalues$p)>1)  # Reduce the number of SNPs for plotting

pdf("Manhattan_plot.pdf",height=6,width=10)
plot(-log10(SNP_reduced$p),pch=16,col=SNP_reduced$color,cex=0.8)
dev.off()

