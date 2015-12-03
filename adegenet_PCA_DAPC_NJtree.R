# Script to explore Aspen Radseq data from Upendra
# Dan Fulop 10/22/2015
# Ideas:
# - lump indiv. into populations by lat/long coords. and also gazing at map
#     - maybe use a certain lat/long radius for lumping localities into same population
#     - map samples in Google Earth and/or ggmap / ggvis + shiny to visualize them? (avoid if possible, *rabbit hole*)

library(ape)
library(adegenet)
library(pegas)
library(seqinr)
library(phangorn)
library(ggplot2)
library(stringr)
library(parallel)
library(genetics)
library(StAMPP)
library(poppr)
library(plyr)


#-------
setwd("~/UCD/Miscelaneous/Upendra_Aspen_demography/")
alle8b <- read.PLINK("/Users/Dani/UCD/Miscelaneous/Upendra_Aspen_demography/HapMap.hmp.8base_filtered.mod.renamed.plk.raw", parallel=T, n.cores=4)
save(alle8b, file="alle8b.Rdata")
load("alle8b.Rdata")

# Calculate a genetic distance matrix of pairwise distances among individuals using Nei's distance

neiD.mat <- stamppNeisD(alle8b, pop=F)
colnames(neiD.mat) <- rownames(neiD.mat) # label columns
save(neiD.mat, file="neiD.mat.Rdata")
write.csv(neiD.mat, file="neiD.mat.csv", quote=F)

# Infer NJ tree with bootstraps
boot.tree <- aboot(alle8b, tree="nj", distance=function(alle8b) as.dist(stamppNeisD(alle8b, pop=F)), sample=100, cutoff=50, showtree=T, quiet=F, root=F)
save(boot.tree, file="boot.tree.Rdata")
write.tree(boot.tree, file="boot_tree_8b.tre")

load("boot.tree.Rdata")
#-----
# pdf(file="allele_data_matrix_8b.pdf")
plot(alle8b) # plot allele matrix -- saved as PNG instead
# dev.off()

# compute and plot allele frequency
myFreq <- glMean(alle8b)
myFreq <- c(myFreq, 1-myFreq)
pdf(file="allele_data_matrix_8b.pdf")
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20, ylim=c(0,3.6))
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)
dev.off()

# pca8b <- glPca(alle8b, parallel=T, n.cores=4, nf=9) # selected 9 axes judging from the eigenvalues plot
pca8b <- glPca(alle8b, parallel=T, n.cores=3) # retained 20 axes
save(pca8b, file="pca8b.Rdata")

# pdf(file="PCA_8b_PC1n2_textLabels.pdf")
# scatter(pca8b, xax=1, yax=2, posi="none", sub="X=PC1, Y=PC2")
# dev.off()
# 
# pdf(file="PCA_8b_PC3n2_textLabels.pdf")
# scatter(pca8b, xax=3, yax=2, posi="none", sub="X=PC3, Y=PC2")
# dev.off()
# 
# pdf(file="PCA_8b_PC1n2_MVcolors.pdf")
myCol <- colorplot(pca8b$scores, pca8b$scores, transp=TRUE, cex=4)
# abline(h=0,v=0, col="grey")
# dev.off()
# 
# pdf(file="PCA_8b_PC3n2_MVcolors.pdf", width=8, height=8)
# colorplot(pca8b$scores[,2:3], pca8b$scores, transp=TRUE, cex=4)
# abline(h=0,v=0, col="grey")
# dev.off()

# add.scatter.eig(pca8b$eig[1:40],2,1,2, posi="topleft", inset=.05, ratio=.3)

load("boot.tree.Rdata")
# plot(boot.tree, type="fan", show.tip=FALSE)
# tiplabels(pch=20, col=myCol, cex=4)
# title("Aspen NJ tree with PCA coloring")

tiplabs <- boot.tree$tip.label
# FIX 2 extra Unknown locations
agr <- grep("_Angelsrest_", tiplabs)
tiplabs[agr] <- str_replace(tiplabs[agr], "_F_", "_Unknown_")

wll <- grep("_WallaWalla_", tiplabs)
tiplabs[wll] <- str_replace(tiplabs[wll], "_F_", "_Unknown_")

flood.idx <- grep("_F_", tiplabs)
nonflood.idx <- which(!grepl("_F_", tiplabs))

tipcols <- vector(mode="character", length=length(tiplabels)) # tip colors for flood and non-flood
tipcols[flood.idx] <- "red"
tipcols[nonflood.idx] <- "black"

pdf(file="NJtree_8b_floodRed.pdf", width=14, height=10)
plot(boot.tree, type="unrooted", no.margin=T, rotate.tree=206, tip.color=tipcols, cex=0.5, lab4ut="axial", font=1)
dev.off()

#-------
load("pca8b.Rdata")
# Give the flood vs non-flood points different shapes
pca.scores <- as.data.frame(pca8b$scores)

# change Angel's Rest and Walla Walla from Flood to Unknown
agr <- grep("_Angelsrest_", rownames(pca.scores))
rownames(pca.scores)[agr] <- str_replace(rownames(pca.scores)[agr], "_F_", "_Unknown_")

wll <- grep("_WallaWalla_", rownames(pca.scores))
rownames(pca.scores)[wll] <- str_replace(rownames(pca.scores)[wll], "_F_", "_Unknown_")

flood <- grep("_F_", rownames(pca.scores))
noflood <- which(!grepl("_F_", rownames(pca.scores)))
unk <- grep("_Unknown_", rownames(pca.scores))
pt.shapes <- vector("character", nrow(pca.scores))
pt.shapes[noflood] <- "1"
pt.shapes[flood] <- "2"
pt.shapes[unk] <- "3"
pca.scores <- cbind(pca.scores, pt.shapes)

pt.shapes <- vector("numeric", nrow(pca.scores))
pt.shapes[noflood] <- 19
pt.shapes[flood] <- 17
pt.shapes[unk] <- 15

# *** Vary dot shape by flood / non-flood ***
pdf(file="NJtree_8b_PCAcoloring.pdf", width=14, height=10)
plot(boot.tree, type="unrooted", show.tip=FALSE, no.margin=T, rotate.tree=206)
tiplabels(col=myCol, cex=2.5, pch=pt.shapes)
dev.off()

# Plot PCs
pc1n2.text <- ggplot(data=pca.scores) + geom_label(aes(x=PC1, y=PC2, label=rownames(pca.scores))) + theme_bw() + coord_equal(xlim=c(-28.5,26))
ggsave(filename="PCA_8b_PC1n2_textLabels.pdf", pc1n2.text, width=14)
pc1n2.col <- ggplot(data=pca.scores, aes(x=PC1, y=PC2)) + geom_point(aes(shape=pt.shapes, size=7), color=myCol) + theme_bw() + coord_equal() + theme(legend.position = "none")
ggsave(filename="PCA_8b_PC1n2_MVcolors.pdf", pc1n2.col, width=14)

pc1n3.text <- ggplot(data=pca.scores) + geom_label(aes(x=PC1, y=PC3, label=rownames(pca.scores))) + theme_bw() + coord_equal(xlim=c(-29.5,25.5))
ggsave(filename="PCA_8b_PC1n3_textLabels.pdf", pc1n3.text, width=14)
pc1n3.col <- ggplot(data=pca.scores, aes(x=PC1, y=PC3)) + geom_point(aes(shape=pt.shapes, size=7), color=myCol) + theme_bw() + coord_equal() + theme(legend.position = "none")
ggsave(filename="PCA_8b_PC1n3_MVcolors.pdf", pc1n3.col, width=14)

pc2n3.text <- ggplot(data=pca.scores) + geom_label(aes(x=PC2, y=PC3, label=rownames(pca.scores))) + theme_bw() + coord_equal(xlim=c(-23,21))
ggsave(filename="PCA_8b_PC2n3_textLabels.pdf", pc2n3.text, height=14)
pc2n3.col <- ggplot(data=pca.scores, aes(x=PC2, y=PC3)) + geom_point(aes(shape=pt.shapes, size=7), color=myCol) + theme_bw() + coord_equal() + theme(legend.position = "none")
ggsave(filename="PCA_8b_PC2n3_MVcolors.pdf", pc2n3.col, height=14)


var.exp <- cumsum(pca8b$eig / sum(pca8b$eig)) # proportion of variance explained by PCs
# PCs 1-3 = ~17% var
# PCs 1-9 = ~25% var

# TRY all 3 criteria ...maybe just choose ~4 clusters
# BTW, using 9 PCs is equivalent to ~25% of the variance
# Since it's based on K-means clustering, I reran the clustering several times until I got a good BIC curve for a given cluster number.
# 7B seems to be the best clustering.
clust.8b.bic <- find.clusters(alle8b, n.pca=9, stat='BIC', max.n.clust=11, glPca=pca8b) # chose 6
str(clust.8b.bic)
clust.8b.bic
names(clust.8b.bic$grp)

agr <- grep("_Angelsrest_", names(clust.8b.bic$grp))
names(clust.8b.bic$grp)[agr] <- str_replace(names(clust.8b.bic$grp)[agr], "_F_", "_Unknown_")

wll <- grep("_WallaWalla_", names(clust.8b.bic$grp))
names(clust.8b.bic$grp)[wll] <- str_replace(names(clust.8b.bic$grp)[wll], "_F_", "_Unknown_")

setdiff(names(clust.8b.bic$grp), rownames(pca.scores)) # FIX the flood classif of Angel's rest and Walla Walla

clust.8b.bic_4clust <- clust.8b.bic
clust.8b.bic_5clust <- clust.8b.bic
clust.8b.bic_6clust <- clust.8b.bic
clust.8b.bic_7clust <- clust.8b.bic
clust.8b.bic_7clustB <- clust.8b.bic

# USE convex hull instead of ellipse!!  ...AND use Structure clusters

pca.scores <- cbind(pca.scores, cl4=clust.8b.bic_4clust$grp, cl5=clust.8b.bic_5clust$grp, cl6=clust.8b.bic_6clust$grp, cl7a=clust.8b.bic_7clust$grp, cl7b=clust.8b.bic_7clustB$grp)

# find.hull <- function(df) df[chull(df$PC1, df$PC2), ] # function that uses chull() to find the convex hull of each cluster
# ddply(pca.scores, "cl7b")
# ddply(pca.scores, "cl7b", find.hull)

# Going old-school w/ split - apply - combine because ddply is failing me
split.pcs <- split(pca.scores[1:3], pca.scores$cl7b)
hulls12.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC1, df$PC2), 1:2]
})
hulls12 <- do.call(rbind, hulls12.list)
hulls12$cluster <- substr(rownames(hulls12), 1, 1)

hulls13.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC1, df$PC3), c(1,3)]
})
hulls13 <- do.call(rbind, hulls13.list)
hulls13$cluster <- substr(rownames(hulls13), 1, 1)

hulls23.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC2, df$PC3), 2:3]
})
hulls23 <- do.call(rbind, hulls23.list)
hulls23$cluster <- substr(rownames(hulls23), 1, 1)

pc1n2.kmClusters7b <- ggplot(data=pca.scores, aes(x=PC1, y=PC2)) + geom_point(aes(shape=pt.shapes, size=7), color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls12[hulls12$cluster=="1",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12[hulls12$cluster=="2",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12[hulls12$cluster=="3",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12[hulls12$cluster=="4",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12[hulls12$cluster=="5",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12[hulls12$cluster=="6",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12[hulls12$cluster=="7",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA)
ggsave(filename="pc1n2.kmClusters7b.pdf", pc1n2.kmClusters7b, width=14)  

pc1n3.kmClusters7b <- ggplot(data=pca.scores, aes(x=PC1, y=PC3)) + geom_point(aes(shape=pt.shapes, size=7), color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls13[hulls13$cluster=="1",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="2",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="3",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="4",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="5",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="6",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="7",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA)
ggsave(filename="pc1n3.kmClusters7b.pdf", pc1n3.kmClusters7b, width=14)  

pc2n3.kmClusters7b <- ggplot(data=pca.scores, aes(x=PC2, y=PC3)) + geom_point(aes(shape=pt.shapes, size=7), color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls23[hulls23$cluster=="1",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23[hulls23$cluster=="2",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23[hulls23$cluster=="3",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23[hulls23$cluster=="4",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23[hulls23$cluster=="5",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23[hulls23$cluster=="6",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23[hulls23$cluster=="7",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA)
ggsave(filename="pc2n3.kmClusters7b.pdf", pc2n3.kmClusters7b, width=14)  

save(pca.scores, file="pca.scores_8b.Rdata")

load("pca.scores_8b.Rdata")

pc.sub <- pca.scores[c(1:3,30)] # sub-matrix of PCA scores that only has the PCs 1-3 and the 7B k-means clustering

# Input fastStructure clusters to display along with PCA
clust.prob <- read.delim("~/UCD/Miscelaneous/Upendra_Aspen_demography/8bp_fastStructure_output/HapMap.hmp.8base_filtered.mod.renamed.output_simple.4.meanQ", 
                         row.names=rownames(pca.scores), header=F)
colnames(clust.prob) <- paste0("c", 1:4)
head(clust.prob)
summary(clust.prob)

clust.prob <- cbind(clust.prob, cl=apply(clust.prob, 1, function(x) which.max(x) )) # Add a cluster assignment

# Cluster uncertainty
summary(apply(clust.prob[1:4], 1, function(x) max(x) < 0.95)) # FALSE:149  TRUE:33
summary(apply(clust.prob[1:4], 1, function(x) max(x) < 0.90)) # FALSE:157  TRUE:25
summary(apply(clust.prob[1:4], 1, function(x) max(x) < 0.80)) # FALSE:166  TRUE:16

clust.prob <- cbind(clust.prob, unc.cl=apply(clust.prob[1:4], 1, function(x) max(x) < 0.80) ) # Add a cluster uncertainty
# save(clust.prob, file="clust.prob.Rdata")

# Combine fastStructure and PCA
str.pc <- cbind(pc.sub, clust.prob)
head(str.pc)

split.pcs <- split(str.pc[1:3], str.pc$cl)
hulls12.fStr.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC1, df$PC2), 1:2]
})
hulls12.fStr <- do.call(rbind, hulls12.fStr.list)
hulls12.fStr$cluster <- substr(rownames(hulls12.fStr), 1, 1)

hulls13.fStr.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC1, df$PC3), c(1,3)]
})
hulls13.fStr <- do.call(rbind, hulls13.fStr.list)
hulls13.fStr$cluster <- substr(rownames(hulls13.fStr), 1, 1)

hulls23.fStr.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC2, df$PC3), 2:3]
})
hulls23.fStr <- do.call(rbind, hulls23.fStr.list)
hulls23.fStr$cluster <- substr(rownames(hulls23.fStr), 1, 1)

# save(str.pc, file="str.pc.Rdata")

# str.pc <- cbind(str.pc, as.character(pt.shapes))

pc1n2.fStr <- ggplot(data=str.pc, aes(x=PC1, y=PC2)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls12.fStr[hulls12.fStr$cluster=="1",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12.fStr[hulls12.fStr$cluster=="2",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12.fStr[hulls12.fStr$cluster=="3",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12.fStr[hulls12.fStr$cluster=="4",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA)
ggsave(filename="pc1n2.fStr.pdf", pc1n2.fStr, width=14)  

pc1n3.fStr <- ggplot(data=str.pc, aes(x=PC1, y=PC3)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="1",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="2",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="3",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="4",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA)
ggsave(filename="pc1n3.fStr.pdf", pc1n3.fStr, width=14)  

pc2n3.fStr <- ggplot(data=str.pc, aes(x=PC2, y=PC3)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls23.fStr[hulls23.fStr$cluster=="1",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23.fStr[hulls23.fStr$cluster=="2",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23.fStr[hulls23.fStr$cluster=="3",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23.fStr[hulls23.fStr$cluster=="4",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA)
ggsave(filename="pc2n3.fStr.pdf", pc2n3.fStr, width=14)  

summary(as.factor(as.character(str.pc$cl)))

# Make a sample-ID label to plot in PC space in order to help with making a curated set of clusters
str.pc$sampleID <- sapply(rownames(str.pc), function(x) str_match(x, "(.+)_.+_.+_.+")[2])

ggplot(data=str.pc, aes(x=PC1, y=PC3)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") + geom_text(aes(label=cl7b), size=3, alpha=0.5) + # + geom_text(aes(label=sampleID), size=3, alpha=0.5) +
  geom_polygon(data=hulls13[hulls13$cluster=="1",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="2",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="3",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="4",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="5",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="6",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13[hulls13$cluster=="7",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA)

ggplot(data=str.pc, aes(x=PC1, y=PC3)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") + geom_text(aes(label=cl), size=3, alpha=0.5) + # + geom_text(aes(label=sampleID), size=3, alpha=0.5) +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="1",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="2",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="3",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.fStr[hulls13.fStr$cluster=="4",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA)

# Create a new modified / curated cluster column
str.pc$mod.cl <- str.pc$cl
# Grab the k-means cluster 6 samples, minus GCB5 and the 2 outlier Kittitas County samples, and assign to a new cluster == 5
str.pc$mod.cl[str.pc$cl7b==6 & str.pc$sampleID!="GCB5" & str.pc$sampleID!="WWA26" & str.pc$sampleID!="WWA03"] <- 5

# Plot the new modified / curated cluster classification
split.pcs <- split(str.pc[1:3], str.pc$mod.cl)
hulls12.modCl.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC1, df$PC2), 1:2]
})
hulls12.modCl <- do.call(rbind, hulls12.modCl.list)
hulls12.modCl$cluster <- substr(rownames(hulls12.modCl), 1, 1)

hulls13.modCl.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC1, df$PC3), c(1,3)]
})
hulls13.modCl <- do.call(rbind, hulls13.modCl.list)
hulls13.modCl$cluster <- substr(rownames(hulls13.modCl), 1, 1)

hulls23.modCl.list <- lapply(split.pcs, function(df) {
  df[chull(df$PC2, df$PC3), 2:3]
})
hulls23.modCl <- do.call(rbind, hulls23.modCl.list)
hulls23.modCl$cluster <- substr(rownames(hulls23.modCl), 1, 1)

pc1n2.modCl <- ggplot(data=str.pc, aes(x=PC1, y=PC2)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls12.modCl[hulls12.modCl$cluster=="1",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12.modCl[hulls12.modCl$cluster=="2",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12.modCl[hulls12.modCl$cluster=="3",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12.modCl[hulls12.modCl$cluster=="4",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls12.modCl[hulls12.modCl$cluster=="5",], aes(x=PC1, y=PC2), color="black", linetype=2, fill=NA) +
  geom_text(x=-16,y=10, hjust="left",label="Willamette + Olympia", size=7) +
  geom_text(x=-16,y=-0.5, hjust="left",label="Glade Creek", size=7) +
  geom_text(x=-7,y=-13, hjust="left",label="Cascade + Sierra", size=7) +
  geom_text(x=3.5,y=-5.5, hjust="left",label="Western Rockies + Eastern WA", size=7) +
  geom_text(x=11,y=3, hjust="right",label="Eastern Rockies +\nCanada + MI", size=7)
ggsave(filename="pc1n2.modCl.pdf", pc1n2.modCl, width=14)  

pc1n3.modCl <- ggplot(data=str.pc, aes(x=PC1, y=PC3)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls13.modCl[hulls13.modCl$cluster=="1",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.modCl[hulls13.modCl$cluster=="2",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.modCl[hulls13.modCl$cluster=="3",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.modCl[hulls13.modCl$cluster=="4",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls13.modCl[hulls13.modCl$cluster=="5",], aes(x=PC1, y=PC3), color="black", linetype=2, fill=NA)
ggsave(filename="pc1n3.modCl.pdf", pc1n3.modCl, width=14)  

pc2n3.modCl <- ggplot(data=str.pc, aes(x=PC2, y=PC3)) + geom_point(aes(size=7), shape=pt.shapes, color=myCol) + theme_bw() + 
  coord_equal() + theme(legend.position = "none") +
  geom_polygon(data=hulls23.modCl[hulls23.modCl$cluster=="1",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23.modCl[hulls23.modCl$cluster=="2",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23.modCl[hulls23.modCl$cluster=="3",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23.modCl[hulls23.modCl$cluster=="4",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA) +
  geom_polygon(data=hulls23.modCl[hulls23.modCl$cluster=="5",], aes(x=PC2, y=PC3), color="black", linetype=2, fill=NA)
ggsave(filename="pc2n3.modCl.pdf", pc2n3.modCl, width=14)  

save(str.pc, file="str.pc_post_modCl.Rdata") # Save the data after the modified cluster creation

#-------------
# READ vignette to see how to display the cluster results
?stat_ellipse

# clust.8b.aic <- find.clusters(alle8b, n.pca=20, stat='AIC', max.n.clust=11, glPca=pca8b) # chose 6
# clust.8b.wss <- find.clusters(alle8b, n.pca=20, stat='WSS', max.n.clust=11, glPca=pca8b) # chose 6

# TO DOs
# 1) Add cluster ellipses to the PCA plots =>
#   1.1) either with stat-ellipse in ggplot and/or w/ ade4.
#   1.2) These can be K-means clusters as well as Structure ones
# 2) Check tree vs. genetic distance with regression
# 3) Do PCoA/ MDS plots??
# 4) Use ggmap / ggvis / shiny to put MVcolored dots on the map
#   4.1) Possibly add cluster polygons as well...
# 5) Do sPCA analysis with the PCA >> multispati workaround? => NO

# IMPORTANT NOTES:
# The Willamette + Olympia cluster could be subdivided into Flood vs. Non-flood, or kept together

# Also, the Calapooia + Angel's rest + Brooks could be carved off as it's own cluster perhaps

# Another good idea may be to remove the Canada+MI samples & also GCB, then do dadi analysis with Rocky, Cascade+Sierra, and Western clusters!!
#  ...with and without western samples divided into Flood / Non-flood