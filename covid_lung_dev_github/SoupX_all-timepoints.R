library(SoupX)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(DropletUtils)
set.seed(33)


#E18_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for E18_c 
E18_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix E18/", gene.column = 2, unique.features = TRUE)
E18_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix E18/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
E18_c.raw.matrix.GE <- E18_c.raw.matrix$`Gene Expression`
E18_c.filtered.matrix.GE <- E18_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
E18_c.sc <- SoupChannel(tod = E18_c.raw.matrix.GE, toc = E18_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for E18_c data
srat_E18_c <-  CreateSeuratObject(E18_c.sc$toc)
srat_E18_c <-  NormalizeData(srat_E18_c)
srat_E18_c <-  FindVariableFeatures(srat_E18_c) 
srat_E18_c <-  ScaleData(srat_E18_c) 
srat_E18_c <-  RunPCA(srat_E18_c,pcs.compute=30) 
srat_E18_c <-  RunUMAP(srat_E18_c,dims=seq(30)) 
srat_E18_c <-  FindNeighbors(srat_E18_c,dims=seq(30)) 
srat_E18_c <-  FindClusters(srat_E18_c,resolution=1) 
E18_c_DR <-  as.data.frame(srat_E18_c@reductions$umap@cell.embeddings)
colnames(E18_c_DR) = c('RD1','RD2')
E18_c_DR$Cluster = factor(srat_E18_c@meta.data[rownames(E18_c_DR),'RNA_snn_res.1'])

#E18_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
E18_c_DR$sftpc = E18_c.sc$toc["Sftpc", rownames(E18_c_DR)]
E18_c_DR$scgb1a1 = E18_c.sc$toc["Scgb1a1", rownames(E18_c_DR)]
E18_c_DR$ptprc = E18_c.sc$toc["Ptprc", rownames(E18_c_DR)]
E18_c_DR$pdgfra = E18_c.sc$toc["Pdgfra", rownames(E18_c_DR)]
E18_c_DR$pecam1 = E18_c.sc$toc["Pecam1", rownames(E18_c_DR)]
E18_c_DR$dcn = E18_c.sc$toc["Dcn", rownames(E18_c_DR)]
ggplot(E18_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(E18_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(E18_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(E18_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(E18_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(E18_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(E18_c.sc, "Scgb1a1", E18_c_DR)
plotMarkerMap(E18_c.sc, "Sftpc", E18_c_DR)
plotMarkerMap(E18_c.sc, "Ptprc", E18_c_DR)
plotMarkerMap(E18_c.sc, "Pdgfra", E18_c_DR)
plotMarkerMap(E18_c.sc, "Pecam1", E18_c_DR)
plotMarkerMap(E18_c.sc, "Dcn", E18_c_DR)

#E18_c, which genes are expressed most highly in the background? 
head(E18_c.sc$soupProfile[order(E18_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(E18_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_E18_c = estimateNonExpressingCells(E18_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(E18_c_DR$Cluster, rownames(E18_c_DR)))
plotMarkerMap(E18_c.sc, geneSet = background_RNA_genes, DR = E18_c_DR, useToEst = useToEst_E18_c)

## calculate contamination and adjust the counts. 
E18_c.sc = setClusters(E18_c.sc, E18_c_DR$Cluster)
E18_c.sc = calculateContaminationFraction(E18_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_E18_c)
head(E18_c.sc$metaData)

## adjust the counts based on contamination fraction
out_E18_c = adjustCounts(E18_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(E18_c.sc$toc > 0)
cntStrained = rowSums(out_E18_c > 0)
mostZeroed_E18_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_E18_c
plotChangeMap(E18_c.sc, out_E18_c, "Cd74", E18_c_DR)
plotChangeMap(E18_c.sc, out_E18_c, "Scgb3a2", E18_c_DR)
plotChangeMap(E18_c.sc, out_E18_c, "Scgb1a1", E18_c_DR)
plotChangeMap(E18_c.sc, out_E18_c, "Sftpc", E18_c_DR)
plotChangeMap(E18_c.sc, out_E18_c, "Dcn", E18_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./E18_c_counts", out_E18_c)

###############################################################


#E18_epi_epi epithelial specific sort SoupX Pipeline

#Read matrix files into R for E18 epithelial cell sorted data
E18_epi_epi.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix epis from EP/", gene.column = 2, unique.features = TRUE)
E18_epi_epi.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix Epis from EP/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
E18_epi_epi.raw.matrix.GE <- E18_epi_epi.raw.matrix$`Gene Expression`
E18_epi_epi.filtered.matrix.GE <- E18_epi_epi.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
E18_epi_epi.sc <- SoupChannel(tod = E18_epi_epi.raw.matrix.GE, toc = E18_epi_epi.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for E18_epi_epi data
srat_E18_epi_epi <-  CreateSeuratObject(E18_epi_epi.sc$toc)
srat_E18_epi_epi <-  NormalizeData(srat_E18_epi_epi)
srat_E18_epi_epi <-  FindVariableFeatures(srat_E18_epi_epi) 
srat_E18_epi_epi <-  ScaleData(srat_E18_epi_epi) 
srat_E18_epi_epi <-  RunPCA(srat_E18_epi_epi,pcs.compute=30) 
srat_E18_epi_epi <-  RunUMAP(srat_E18_epi_epi,dims=seq(30)) 
srat_E18_epi_epi <-  FindNeighbors(srat_E18_epi_epi,dims=seq(30)) 
srat_E18_epi_epi <-  FindClusters(srat_E18_epi_epi,resolution=1) 
E18_epi_epi_DR <-  as.data.frame(srat_E18_epi_epi@reductions$umap@cell.embeddings)
colnames(E18_epi_epi_DR) = c('RD1','RD2')
E18_epi_epi_DR$Cluster = factor(srat_E18_epi_epi@meta.data[rownames(E18_epi_epi_DR),'RNA_snn_res.1'])

# E18_epi_epi Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
E18_epi_epi_DR$sftpc = E18_epi_epi.sc$toc["Sftpc", rownames(E18_epi_epi_DR)]
E18_epi_epi_DR$scgb1a1 = E18_epi_epi.sc$toc["Scgb1a1", rownames(E18_epi_epi_DR)]
E18_epi_epi_DR$ptprc = E18_epi_epi.sc$toc["Ptprc", rownames(E18_epi_epi_DR)]
E18_epi_epi_DR$pdgfra = E18_epi_epi.sc$toc["Pdgfra", rownames(E18_epi_epi_DR)]
E18_epi_epi_DR$pecam1 = E18_epi_epi.sc$toc["Pecam1", rownames(E18_epi_epi_DR)]
E18_epi_epi_DR$dcn = E18_epi_epi.sc$toc["Dcn", rownames(E18_epi_epi_DR)]
ggplot(E18_epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(E18_epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(E18_epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(E18_epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(E18_epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(E18_epi_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(E18_epi_epi.sc, "Scgb1a1", E18_epi_epi_DR)
plotMarkerMap(E18_epi_epi.sc, "Sftpc", E18_epi_epi_DR)
plotMarkerMap(E18_epi_epi.sc, "Ptprc", E18_epi_epi_DR)
plotMarkerMap(E18_epi_epi.sc, "Pdgfra", E18_epi_epi_DR)
plotMarkerMap(E18_epi_epi.sc, "Pecam1", E18_epi_epi_DR)
plotMarkerMap(E18_epi_epi.sc, "Dcn", E18_epi_epi_DR)

# E18_epi_epi, which genes are expressed most highly in the background? 
head(E18_epi_epi.sc$soupProfile[order(E18_epi_epi.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(E18_epi_epi.sc)

#Specify background RNA genes for non-epithelial cell specific sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_E18_epi_epi = estimateNonExpressingCells(E18_epi_epi.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(E18_epi_epi_DR$Cluster, rownames(E18_epi_epi_DR)))
plotMarkerMap(E18_epi_epi.sc, geneSet = background_RNA_genes, DR = E18_epi_epi_DR, useToEst = useToEst_E18_epi_epi)

## calculate contamination and adjust the counts. 
E18_epi_epi.sc = setClusters(E18_epi_epi.sc, E18_epi_epi_DR$Cluster)
E18_epi_epi.sc = calculateContaminationFraction(E18_epi_epi.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_E18_epi_epi)
head(E18_epi_epi.sc$metaData)

## adjust the counts based on contamination fraction
out_E18_epi_epi = adjustCounts(E18_epi_epi.sc)

## rechecking some plots and genes to check the genes are mostly specific to clusters. 
cntSoggy = rowSums(E18_epi_epi.sc$toc > 0)
cntStrained = rowSums(out_E18_epi_epi > 0)
mostZeroed_E18_epi_epi = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_E18_epi_epi
plotChangeMap(E18_epi_epi.sc, out_E18_epi_epi, "Cd74", E18_epi_epi_DR)
plotChangeMap(E18_epi_epi.sc, out_E18_epi_epi, "Scgb3a2", E18_epi_epi_DR)
plotChangeMap(E18_epi_epi.sc, out_E18_epi_epi, "Scgb1a1", E18_epi_epi_DR)
plotChangeMap(E18_epi_epi.sc, out_E18_epi_epi, "Sftpc", E18_epi_epi_DR)
plotChangeMap(E18_epi_epi.sc, out_E18_epi_epi, "Dcn", E18_epi_epi_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./E18_epi_epi_counts", out_E18_epi_epi)

##########################################################

#P0_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for P0 
P0_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix JS6/", gene.column = 2, unique.features = TRUE)
P0_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix JS6/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
P0_c.raw.matrix.GE <- P0_c.raw.matrix$`Gene Expression`
P0_c.filtered.matrix.GE <- P0_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
P0_c.sc <- SoupChannel(tod = P0_c.raw.matrix.GE, toc = P0_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for P0_c data
srat_P0_c <-  CreateSeuratObject(P0_c.sc$toc)
srat_P0_c <-  NormalizeData(srat_P0_c)
srat_P0_c <-  FindVariableFeatures(srat_P0_c) 
srat_P0_c <-  ScaleData(srat_P0_c) 
srat_P0_c <-  RunPCA(srat_P0_c,pcs.compute=30) 
srat_P0_c <-  RunUMAP(srat_P0_c,dims=seq(30)) 
srat_P0_c <-  FindNeighbors(srat_P0_c,dims=seq(30)) 
srat_P0_c <-  FindClusters(srat_P0_c,resolution=1) 
P0_c_DR <-  as.data.frame(srat_P0_c@reductions$umap@cell.embeddings)
colnames(P0_c_DR) = c('RD1','RD2')
P0_c_DR$Cluster = factor(srat_P0_c@meta.data[rownames(P0_c_DR),'RNA_snn_res.1'])

#P0_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
P0_c_DR$sftpc = P0_c.sc$toc["Sftpc", rownames(P0_c_DR)]
P0_c_DR$scgb1a1 = P0_c.sc$toc["Scgb1a1", rownames(P0_c_DR)]
P0_c_DR$ptprc = P0_c.sc$toc["Ptprc", rownames(P0_c_DR)]
P0_c_DR$pdgfra = P0_c.sc$toc["Pdgfra", rownames(P0_c_DR)]
P0_c_DR$pecam1 = P0_c.sc$toc["Pecam1", rownames(P0_c_DR)]
P0_c_DR$dcn = P0_c.sc$toc["Dcn", rownames(P0_c_DR)]
ggplot(P0_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(P0_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(P0_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(P0_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(P0_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(P0_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(P0_c.sc, "Scgb1a1", P0_c_DR)
plotMarkerMap(P0_c.sc, "Sftpc", P0_c_DR)
plotMarkerMap(P0_c.sc, "Ptprc", P0_c_DR)
plotMarkerMap(P0_c.sc, "Pdgfra", P0_c_DR)
plotMarkerMap(P0_c.sc, "Pecam1", P0_c_DR)
plotMarkerMap(P0_c.sc, "Dcn", P0_c_DR)

#P0_c, which genes are expressed most highly in the background? 
head(P0_c.sc$soupProfile[order(P0_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(P0_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_P0_c = estimateNonExpressingCells(P0_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(P0_c_DR$Cluster, rownames(P0_c_DR)))
plotMarkerMap(P0_c.sc, geneSet = background_RNA_genes, DR = P0_c_DR, useToEst = useToEst_P0_c)

## calculate contamination and adjust the counts. 
P0_c.sc = setClusters(P0_c.sc, P0_c_DR$Cluster)
P0_c.sc = calculateContaminationFraction(P0_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_P0_c)
head(P0_c.sc$metaData)

## adjust the counts based on contamination fraction
out_P0_c = adjustCounts(P0_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(P0_c.sc$toc > 0)
cntStrained = rowSums(out_P0_c > 0)
mostZeroed_P0_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_P0_c
plotChangeMap(P0_c.sc, out_P0_c, "Cd74", P0_c_DR)
plotChangeMap(P0_c.sc, out_P0_c, "Scgb3a2", P0_c_DR)
plotChangeMap(P0_c.sc, out_P0_c, "Scgb1a1", P0_c_DR)
plotChangeMap(P0_c.sc, out_P0_c, "Sftpc", P0_c_DR)
plotChangeMap(P0_c.sc, out_P0_c, "Dcn", P0_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./P0_c_counts", out_P0_c)

###############################################################

#P0_2_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for P0_2_c
P0_2_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix JS-1 P0/", gene.column = 2, unique.features = TRUE)
P0_2_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix JS-1 P0/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
P0_2_c.raw.matrix.GE <- P0_2_c.raw.matrix$`Gene Expression`
P0_2_c.filtered.matrix.GE <- P0_2_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
P0_2_c.sc <- SoupChannel(tod = P0_2_c.raw.matrix.GE, toc = P0_2_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for P0_2_c data
srat_P0_2_c <-  CreateSeuratObject(P0_2_c.sc$toc)
srat_P0_2_c <-  NormalizeData(srat_P0_2_c)
srat_P0_2_c <-  FindVariableFeatures(srat_P0_2_c) 
srat_P0_2_c <-  ScaleData(srat_P0_2_c) 
srat_P0_2_c <-  RunPCA(srat_P0_2_c,pcs.compute=30) 
srat_P0_2_c <-  RunUMAP(srat_P0_2_c,dims=seq(30)) 
srat_P0_2_c <-  FindNeighbors(srat_P0_2_c,dims=seq(30)) 
srat_P0_2_c <-  FindClusters(srat_P0_2_c,resolution=1) 
P0_2_c_DR <-  as.data.frame(srat_P0_2_c@reductions$umap@cell.embeddings)
colnames(P0_2_c_DR) = c('RD1','RD2')
P0_2_c_DR$Cluster = factor(srat_P0_2_c@meta.data[rownames(P0_2_c_DR),'RNA_snn_res.1'])

#P0_2_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
P0_2_c_DR$sftpc = P0_2_c.sc$toc["Sftpc", rownames(P0_2_c_DR)]
P0_2_c_DR$scgb1a1 = P0_2_c.sc$toc["Scgb1a1", rownames(P0_2_c_DR)]
P0_2_c_DR$ptprc = P0_2_c.sc$toc["Ptprc", rownames(P0_2_c_DR)]
P0_2_c_DR$pdgfra = P0_2_c.sc$toc["Pdgfra", rownames(P0_2_c_DR)]
P0_2_c_DR$pecam1 = P0_2_c.sc$toc["Pecam1", rownames(P0_2_c_DR)]
P0_2_c_DR$dcn = P0_2_c.sc$toc["Dcn", rownames(P0_2_c_DR)]
ggplot(P0_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(P0_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(P0_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(P0_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(P0_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(P0_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(P0_2_c.sc, "Scgb1a1", P0_2_c_DR)
plotMarkerMap(P0_2_c.sc, "Sftpc", P0_2_c_DR)
plotMarkerMap(P0_2_c.sc, "Ptprc", P0_2_c_DR)
plotMarkerMap(P0_2_c.sc, "Pdgfra", P0_2_c_DR)
plotMarkerMap(P0_2_c.sc, "Pecam1", P0_2_c_DR)
plotMarkerMap(P0_2_c.sc, "Dcn", P0_2_c_DR)

#P0_2_c, which genes are expressed most highly in the background? 
head(P0_2_c.sc$soupProfile[order(P0_2_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(P0_2_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_P0_2_c = estimateNonExpressingCells(P0_2_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(P0_2_c_DR$Cluster, rownames(P0_2_c_DR)))
plotMarkerMap(P0_2_c.sc, geneSet = background_RNA_genes, DR = P0_2_c_DR, useToEst = useToEst_P0_2_c)

## calculate contamination and adjust the counts. 
P0_2_c.sc = setClusters(P0_2_c.sc, P0_2_c_DR$Cluster)
P0_2_c.sc = calculateContaminationFraction(P0_2_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_P0_2_c)
head(P0_2_c.sc$metaData)

## adjust the counts based on contamination fraction
out_P0_2_c = adjustCounts(P0_2_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(P0_2_c.sc$toc > 0)
cntStrained = rowSums(out_P0_2_c > 0)
mostZeroed_P0_2_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_P0_2_c
plotChangeMap(P0_2_c.sc, out_P0_2_c, "Cd74", P0_2_c_DR)
plotChangeMap(P0_2_c.sc, out_P0_2_c, "Scgb3a2", P0_2_c_DR)
plotChangeMap(P0_2_c.sc, out_P0_2_c, "Scgb1a1", P0_2_c_DR)
plotChangeMap(P0_2_c.sc, out_P0_2_c, "Sftpc", P0_2_c_DR)
plotChangeMap(P0_2_c.sc, out_P0_2_c, "Dcn", P0_2_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./P0_2_c_counts", out_P0_2_c)

###############################################################

#P7_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for P7_c
P7_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix P7/", gene.column = 2, unique.features = TRUE)
P7_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix P7/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
P7_c.raw.matrix.GE <- P7_c.raw.matrix$`Gene Expression`
P7_c.filtered.matrix.GE <- P7_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
P7_c.sc <- SoupChannel(tod = P7_c.raw.matrix.GE, toc = P7_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for P7_c data
srat_P7_c <-  CreateSeuratObject(P7_c.sc$toc)
srat_P7_c <-  NormalizeData(srat_P7_c)
srat_P7_c <-  FindVariableFeatures(srat_P7_c) 
srat_P7_c <-  ScaleData(srat_P7_c) 
srat_P7_c <-  RunPCA(srat_P7_c,pcs.compute=30) 
srat_P7_c <-  RunUMAP(srat_P7_c,dims=seq(30)) 
srat_P7_c <-  FindNeighbors(srat_P7_c,dims=seq(30)) 
srat_P7_c <-  FindClusters(srat_P7_c,resolution=1) 
P7_c_DR <-  as.data.frame(srat_P7_c@reductions$umap@cell.embeddings)
colnames(P7_c_DR) = c('RD1','RD2')
P7_c_DR$Cluster = factor(srat_P7_c@meta.data[rownames(P7_c_DR),'RNA_snn_res.1'])

#P7_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
P7_c_DR$sftpc = P7_c.sc$toc["Sftpc", rownames(P7_c_DR)]
P7_c_DR$scgb1a1 = P7_c.sc$toc["Scgb1a1", rownames(P7_c_DR)]
P7_c_DR$ptprc = P7_c.sc$toc["Ptprc", rownames(P7_c_DR)]
P7_c_DR$pdgfra = P7_c.sc$toc["Pdgfra", rownames(P7_c_DR)]
P7_c_DR$pecam1 = P7_c.sc$toc["Pecam1", rownames(P7_c_DR)]
P7_c_DR$dcn = P7_c.sc$toc["Dcn", rownames(P7_c_DR)]
ggplot(P7_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(P7_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(P7_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(P7_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(P7_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(P7_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(P7_c.sc, "Scgb1a1", P7_c_DR)
plotMarkerMap(P7_c.sc, "Sftpc", P7_c_DR)
plotMarkerMap(P7_c.sc, "Ptprc", P7_c_DR)
plotMarkerMap(P7_c.sc, "Pdgfra", P7_c_DR)
plotMarkerMap(P7_c.sc, "Pecam1", P7_c_DR)
plotMarkerMap(P7_c.sc, "Dcn", P7_c_DR)

#P7_c, which genes are expressed most highly in the background? 
head(P7_c.sc$soupProfile[order(P7_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(P7_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_P7_c = estimateNonExpressingCells(P7_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(P7_c_DR$Cluster, rownames(P7_c_DR)))
plotMarkerMap(P7_c.sc, geneSet = background_RNA_genes, DR = P7_c_DR, useToEst = useToEst_P7_c)

## calculate contamination and adjust the counts. 
P7_c.sc = setClusters(P7_c.sc, P7_c_DR$Cluster)
P7_c.sc = calculateContaminationFraction(P7_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_P7_c)
head(P7_c.sc$metaData)

## adjust the counts based on contamination fraction
out_P7_c = adjustCounts(P7_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(P7_c.sc$toc > 0)
cntStrained = rowSums(out_P7_c > 0)
mostZeroed_P7_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_P7_c
plotChangeMap(P7_c.sc, out_P7_c, "Cd74", P7_c_DR)
plotChangeMap(P7_c.sc, out_P7_c, "Scgb3a2", P7_c_DR)
plotChangeMap(P7_c.sc, out_P7_c, "Scgb1a1", P7_c_DR)
plotChangeMap(P7_c.sc, out_P7_c, "Sftpc", P7_c_DR)
plotChangeMap(P7_c.sc, out_P7_c, "Dcn", P7_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./P7_c_counts", out_P7_c)

###############################################################

#P7_3_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for P7_3_c
P7_3_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix P7-3 JS13/", gene.column = 2, unique.features = TRUE)
P7_3_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix P7-3 JS13/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
P7_3_c.raw.matrix.GE <- P7_3_c.raw.matrix$`Gene Expression`
P7_3_c.filtered.matrix.GE <- P7_3_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
P7_3_c.sc <- SoupChannel(tod = P7_3_c.raw.matrix.GE, toc = P7_3_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for P7_3_c data
srat_P7_3_c <-  CreateSeuratObject(P7_3_c.sc$toc)
srat_P7_3_c <-  NormalizeData(srat_P7_3_c)
srat_P7_3_c <-  FindVariableFeatures(srat_P7_3_c) 
srat_P7_3_c <-  ScaleData(srat_P7_3_c) 
srat_P7_3_c <-  RunPCA(srat_P7_3_c,pcs.compute=30) 
srat_P7_3_c <-  RunUMAP(srat_P7_3_c,dims=seq(30)) 
srat_P7_3_c <-  FindNeighbors(srat_P7_3_c,dims=seq(30)) 
srat_P7_3_c <-  FindClusters(srat_P7_3_c,resolution=1) 
P7_3_c_DR <-  as.data.frame(srat_P7_3_c@reductions$umap@cell.embeddings)
colnames(P7_3_c_DR) = c('RD1','RD2')
P7_3_c_DR$Cluster = factor(srat_P7_3_c@meta.data[rownames(P7_3_c_DR),'RNA_snn_res.1'])

#P7_3_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
P7_3_c_DR$sftpc = P7_3_c.sc$toc["Sftpc", rownames(P7_3_c_DR)]
P7_3_c_DR$scgb1a1 = P7_3_c.sc$toc["Scgb1a1", rownames(P7_3_c_DR)]
P7_3_c_DR$ptprc = P7_3_c.sc$toc["Ptprc", rownames(P7_3_c_DR)]
P7_3_c_DR$pdgfra = P7_3_c.sc$toc["Pdgfra", rownames(P7_3_c_DR)]
P7_3_c_DR$pecam1 = P7_3_c.sc$toc["Pecam1", rownames(P7_3_c_DR)]
P7_3_c_DR$dcn = P7_3_c.sc$toc["Dcn", rownames(P7_3_c_DR)]
ggplot(P7_3_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(P7_3_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(P7_3_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(P7_3_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(P7_3_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(P7_3_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(P7_3_c.sc, "Scgb1a1", P7_3_c_DR)
plotMarkerMap(P7_3_c.sc, "Sftpc", P7_3_c_DR)
plotMarkerMap(P7_3_c.sc, "Ptprc", P7_3_c_DR)
plotMarkerMap(P7_3_c.sc, "Pdgfra", P7_3_c_DR)
plotMarkerMap(P7_3_c.sc, "Pecam1", P7_3_c_DR)
plotMarkerMap(P7_3_c.sc, "Dcn", P7_3_c_DR)

#P7_3_c, which genes are expressed most highly in the background? 
head(P7_3_c.sc$soupProfile[order(P7_3_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(P7_3_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_P7_3_c = estimateNonExpressingCells(P7_3_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(P7_3_c_DR$Cluster, rownames(P7_3_c_DR)))
plotMarkerMap(P7_3_c.sc, geneSet = background_RNA_genes, DR = P7_3_c_DR, useToEst = useToEst_P7_3_c)

## calculate contamination and adjust the counts. 
P7_3_c.sc = setClusters(P7_3_c.sc, P7_3_c_DR$Cluster)
P7_3_c.sc = calculateContaminationFraction(P7_3_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_P7_3_c)
head(P7_3_c.sc$metaData)

## adjust the counts based on contamination fraction
out_P7_3_c = adjustCounts(P7_3_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(P7_3_c.sc$toc > 0)
cntStrained = rowSums(out_P7_3_c > 0)
mostZeroed_P7_3_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_P7_3_c
plotChangeMap(P7_3_c.sc, out_P7_3_c, "Cd74", P7_3_c_DR)
plotChangeMap(P7_3_c.sc, out_P7_3_c, "Scgb3a2", P7_3_c_DR)
plotChangeMap(P7_3_c.sc, out_P7_3_c, "Scgb1a1", P7_3_c_DR)
plotChangeMap(P7_3_c.sc, out_P7_3_c, "Sftpc", P7_3_c_DR)
plotChangeMap(P7_3_c.sc, out_P7_3_c, "Dcn", P7_3_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./P7_3_c_counts", out_P7_3_c)

###############################################################

#P7_enriched epithelial specific SoupX Pipeline

#Read matrix files into R for P7_enriched epithelial cell sorted data
P7_enriched_epi.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix P7 enriched/", gene.column = 2, unique.features = TRUE)
P7_enriched_epi.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix P7 enriched/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
P7_enriched_epi.raw.matrix.GE <- P7_enriched_epi.raw.matrix$`Gene Expression`
P7_enriched_epi.filtered.matrix.GE <- P7_enriched_epi.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
P7_enriched_epi.sc <- SoupChannel(tod = P7_enriched_epi.raw.matrix.GE, toc = P7_enriched_epi.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for P7_enriched_epi data
srat_P7_enriched_epi <-  CreateSeuratObject(P7_enriched_epi.sc$toc)
srat_P7_enriched_epi <-  NormalizeData(srat_P7_enriched_epi)
srat_P7_enriched_epi <-  FindVariableFeatures(srat_P7_enriched_epi) 
srat_P7_enriched_epi <-  ScaleData(srat_P7_enriched_epi) 
srat_P7_enriched_epi <-  RunPCA(srat_P7_enriched_epi,pcs.compute=30) 
srat_P7_enriched_epi <-  RunUMAP(srat_P7_enriched_epi,dims=seq(30)) 
srat_P7_enriched_epi <-  FindNeighbors(srat_P7_enriched_epi,dims=seq(30)) 
srat_P7_enriched_epi <-  FindClusters(srat_P7_enriched_epi,resolution=1) 
P7_enriched_epi_DR <-  as.data.frame(srat_P7_enriched_epi@reductions$umap@cell.embeddings)
colnames(P7_enriched_epi_DR) = c('RD1','RD2')
P7_enriched_epi_DR$Cluster = factor(srat_P7_enriched_epi@meta.data[rownames(P7_enriched_epi_DR),'RNA_snn_res.1'])

# P7_enriched_epi Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
P7_enriched_epi_DR$sftpc = P7_enriched_epi.sc$toc["Sftpc", rownames(P7_enriched_epi_DR)]
P7_enriched_epi_DR$scgb1a1 = P7_enriched_epi.sc$toc["Scgb1a1", rownames(P7_enriched_epi_DR)]
P7_enriched_epi_DR$ptprc = P7_enriched_epi.sc$toc["Ptprc", rownames(P7_enriched_epi_DR)]
P7_enriched_epi_DR$pdgfra = P7_enriched_epi.sc$toc["Pdgfra", rownames(P7_enriched_epi_DR)]
P7_enriched_epi_DR$pecam1 = P7_enriched_epi.sc$toc["Pecam1", rownames(P7_enriched_epi_DR)]
P7_enriched_epi_DR$dcn = P7_enriched_epi.sc$toc["Dcn", rownames(P7_enriched_epi_DR)]
ggplot(P7_enriched_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(P7_enriched_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(P7_enriched_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(P7_enriched_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(P7_enriched_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(P7_enriched_epi_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(P7_enriched_epi.sc, "Scgb1a1", P7_enriched_epi_DR)
plotMarkerMap(P7_enriched_epi.sc, "Sftpc", P7_enriched_epi_DR)
plotMarkerMap(P7_enriched_epi.sc, "Ptprc", P7_enriched_epi_DR)
plotMarkerMap(P7_enriched_epi.sc, "Pdgfra", P7_enriched_epi_DR)
plotMarkerMap(P7_enriched_epi.sc, "Pecam1", P7_enriched_epi_DR)
plotMarkerMap(P7_enriched_epi.sc, "Dcn", P7_enriched_epi_DR)

# P7_enriched_epi, which genes are expressed most highly in the background? 
head(P7_enriched_epi.sc$soupProfile[order(P7_enriched_epi.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(P7_enriched_epi.sc)

#Specify background RNA genes for non-epithelial cell specific sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_P7_enriched_epi = estimateNonExpressingCells(P7_enriched_epi.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(P7_enriched_epi_DR$Cluster, rownames(P7_enriched_epi_DR)))
plotMarkerMap(P7_enriched_epi.sc, geneSet = background_RNA_genes, DR = P7_enriched_epi_DR, useToEst = useToEst_P7_enriched_epi)

## calculate contamination and adjust the counts. 
P7_enriched_epi.sc = setClusters(P7_enriched_epi.sc, P7_enriched_epi_DR$Cluster)
P7_enriched_epi.sc = calculateContaminationFraction(P7_enriched_epi.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_P7_enriched_epi)
head(P7_enriched_epi.sc$metaData)

## adjust the counts based on contamination fraction
out_P7_enriched_epi = adjustCounts(P7_enriched_epi.sc)

## rechecking some plots and genes to check the genes are mostly specific to clusters. 
cntSoggy = rowSums(P7_enriched_epi.sc$toc > 0)
cntStrained = rowSums(out_P7_enriched_epi > 0)
mostZeroed_P7_enriched_epi = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_P7_enriched_epi
plotChangeMap(P7_enriched_epi.sc, out_P7_enriched_epi, "Cd74", P7_enriched_epi_DR)
plotChangeMap(P7_enriched_epi.sc, out_P7_enriched_epi, "Scgb3a2", P7_enriched_epi_DR)
plotChangeMap(P7_enriched_epi.sc, out_P7_enriched_epi, "Scgb1a1", P7_enriched_epi_DR)
plotChangeMap(P7_enriched_epi.sc, out_P7_enriched_epi, "Sftpc", P7_enriched_epi_DR)
plotChangeMap(P7_enriched_epi.sc, out_P7_enriched_epi, "Dcn", P7_enriched_epi_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./P7_enriched_epi_counts", out_P7_enriched_epi)

##########################################################

#P14_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for P14_c
P14_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix P14/", gene.column = 2, unique.features = TRUE)
P14_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix P14/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
P14_c.raw.matrix.GE <- P14_c.raw.matrix$`Gene Expression`
P14_c.filtered.matrix.GE <- P14_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
P14_c.sc <- SoupChannel(tod = P14_c.raw.matrix.GE, toc = P14_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for P14_c data
srat_P14_c <-  CreateSeuratObject(P14_c.sc$toc)
srat_P14_c <-  NormalizeData(srat_P14_c)
srat_P14_c <-  FindVariableFeatures(srat_P14_c) 
srat_P14_c <-  ScaleData(srat_P14_c) 
srat_P14_c <-  RunPCA(srat_P14_c,pcs.compute=30) 
srat_P14_c <-  RunUMAP(srat_P14_c,dims=seq(30)) 
srat_P14_c <-  FindNeighbors(srat_P14_c,dims=seq(30)) 
srat_P14_c <-  FindClusters(srat_P14_c,resolution=1) 
P14_c_DR <-  as.data.frame(srat_P14_c@reductions$umap@cell.embeddings)
colnames(P14_c_DR) = c('RD1','RD2')
P14_c_DR$Cluster = factor(srat_P14_c@meta.data[rownames(P14_c_DR),'RNA_snn_res.1'])

#P14_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
P14_c_DR$sftpc = P14_c.sc$toc["Sftpc", rownames(P14_c_DR)]
P14_c_DR$scgb1a1 = P14_c.sc$toc["Scgb1a1", rownames(P14_c_DR)]
P14_c_DR$ptprc = P14_c.sc$toc["Ptprc", rownames(P14_c_DR)]
P14_c_DR$pdgfra = P14_c.sc$toc["Pdgfra", rownames(P14_c_DR)]
P14_c_DR$pecam1 = P14_c.sc$toc["Pecam1", rownames(P14_c_DR)]
P14_c_DR$dcn = P14_c.sc$toc["Dcn", rownames(P14_c_DR)]
ggplot(P14_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(P14_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(P14_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(P14_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(P14_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(P14_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(P14_c.sc, "Scgb1a1", P14_c_DR)
plotMarkerMap(P14_c.sc, "Sftpc", P14_c_DR)
plotMarkerMap(P14_c.sc, "Ptprc", P14_c_DR)
plotMarkerMap(P14_c.sc, "Pdgfra", P14_c_DR)
plotMarkerMap(P14_c.sc, "Pecam1", P14_c_DR)
plotMarkerMap(P14_c.sc, "Dcn", P14_c_DR)

#P14_c, which genes are expressed most highly in the background? 
head(P14_c.sc$soupProfile[order(P14_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(P14_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_P14_c = estimateNonExpressingCells(P14_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(P14_c_DR$Cluster, rownames(P14_c_DR)))
plotMarkerMap(P14_c.sc, geneSet = background_RNA_genes, DR = P14_c_DR, useToEst = useToEst_P14_c)

## calculate contamination and adjust the counts. 
P14_c.sc = setClusters(P14_c.sc, P14_c_DR$Cluster)
P14_c.sc = calculateContaminationFraction(P14_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_P14_c)
head(P14_c.sc$metaData)

## adjust the counts based on contamination fraction
out_P14_c = adjustCounts(P14_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(P14_c.sc$toc > 0)
cntStrained = rowSums(out_P14_c > 0)
mostZeroed_P14_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_P14_c
plotChangeMap(P14_c.sc, out_P14_c, "Cd74", P14_c_DR)
plotChangeMap(P14_c.sc, out_P14_c, "Scgb3a2", P14_c_DR)
plotChangeMap(P14_c.sc, out_P14_c, "Scgb1a1", P14_c_DR)
plotChangeMap(P14_c.sc, out_P14_c, "Sftpc", P14_c_DR)
plotChangeMap(P14_c.sc, out_P14_c, "Dcn", P14_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./P14_c_counts", out_P14_c)

###############################################################

#P14_2_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for P14_2_c
P14_2_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix JS10 P14-2/", gene.column = 2, unique.features = TRUE)
P14_2_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix JS10 P14-2/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
P14_2_c.raw.matrix.GE <- P14_2_c.raw.matrix$`Gene Expression`
P14_2_c.filtered.matrix.GE <- P14_2_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
P14_2_c.sc <- SoupChannel(tod = P14_2_c.raw.matrix.GE, toc = P14_2_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for P14_2_c data
srat_P14_2_c <-  CreateSeuratObject(P14_2_c.sc$toc)
srat_P14_2_c <-  NormalizeData(srat_P14_2_c)
srat_P14_2_c <-  FindVariableFeatures(srat_P14_2_c) 
srat_P14_2_c <-  ScaleData(srat_P14_2_c) 
srat_P14_2_c <-  RunPCA(srat_P14_2_c,pcs.compute=30) 
srat_P14_2_c <-  RunUMAP(srat_P14_2_c,dims=seq(30)) 
srat_P14_2_c <-  FindNeighbors(srat_P14_2_c,dims=seq(30)) 
srat_P14_2_c <-  FindClusters(srat_P14_2_c,resolution=1) 
P14_2_c_DR <-  as.data.frame(srat_P14_2_c@reductions$umap@cell.embeddings)
colnames(P14_2_c_DR) = c('RD1','RD2')
P14_2_c_DR$Cluster = factor(srat_P14_2_c@meta.data[rownames(P14_2_c_DR),'RNA_snn_res.1'])

#P14_2_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
P14_2_c_DR$sftpc = P14_2_c.sc$toc["Sftpc", rownames(P14_2_c_DR)]
P14_2_c_DR$scgb1a1 = P14_2_c.sc$toc["Scgb1a1", rownames(P14_2_c_DR)]
P14_2_c_DR$ptprc = P14_2_c.sc$toc["Ptprc", rownames(P14_2_c_DR)]
P14_2_c_DR$pdgfra = P14_2_c.sc$toc["Pdgfra", rownames(P14_2_c_DR)]
P14_2_c_DR$pecam1 = P14_2_c.sc$toc["Pecam1", rownames(P14_2_c_DR)]
P14_2_c_DR$dcn = P14_2_c.sc$toc["Dcn", rownames(P14_2_c_DR)]
ggplot(P14_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(P14_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(P14_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(P14_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(P14_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(P14_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(P14_2_c.sc, "Scgb1a1", P14_2_c_DR)
plotMarkerMap(P14_2_c.sc, "Sftpc", P14_2_c_DR)
plotMarkerMap(P14_2_c.sc, "Ptprc", P14_2_c_DR)
plotMarkerMap(P14_2_c.sc, "Pdgfra", P14_2_c_DR)
plotMarkerMap(P14_2_c.sc, "Pecam1", P14_2_c_DR)
plotMarkerMap(P14_2_c.sc, "Dcn", P14_2_c_DR)

#P14_2_c, which genes are expressed most highly in the background? 
head(P14_2_c.sc$soupProfile[order(P14_2_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(P14_2_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_P14_2_c = estimateNonExpressingCells(P14_2_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(P14_2_c_DR$Cluster, rownames(P14_2_c_DR)))
plotMarkerMap(P14_2_c.sc, geneSet = background_RNA_genes, DR = P14_2_c_DR, useToEst = useToEst_P14_2_c)

## calculate contamination and adjust the counts. 
P14_2_c.sc = setClusters(P14_2_c.sc, P14_2_c_DR$Cluster)
P14_2_c.sc = calculateContaminationFraction(P14_2_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_P14_2_c)
head(P14_2_c.sc$metaData)

## adjust the counts based on contamination fraction
out_P14_2_c = adjustCounts(P14_2_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(P14_2_c.sc$toc > 0)
cntStrained = rowSums(out_P14_2_c > 0)
mostZeroed_P14_2_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_P14_2_c
plotChangeMap(P14_2_c.sc, out_P14_2_c, "Cd74", P14_2_c_DR)
plotChangeMap(P14_2_c.sc, out_P14_2_c, "Scgb3a2", P14_2_c_DR)
plotChangeMap(P14_2_c.sc, out_P14_2_c, "Scgb1a1", P14_2_c_DR)
plotChangeMap(P14_2_c.sc, out_P14_2_c, "Sftpc", P14_2_c_DR)
plotChangeMap(P14_2_c.sc, out_P14_2_c, "Dcn", P14_2_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./P14_2_c_counts", out_P14_2_c)

###############################################################

#Adult_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for Adult_c
Adult_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix JS5/", gene.column = 2, unique.features = TRUE)
Adult_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix JS5/", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
Adult_c.raw.matrix.GE <- Adult_c.raw.matrix$`Gene Expression`
Adult_c.filtered.matrix.GE <- Adult_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
Adult_c.sc <- SoupChannel(tod = Adult_c.raw.matrix.GE, toc = Adult_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for Adult_c data
srat_Adult_c <-  CreateSeuratObject(Adult_c.sc$toc)
srat_Adult_c <-  NormalizeData(srat_Adult_c)
srat_Adult_c <-  FindVariableFeatures(srat_Adult_c) 
srat_Adult_c <-  ScaleData(srat_Adult_c) 
srat_Adult_c <-  RunPCA(srat_Adult_c,pcs.compute=30) 
srat_Adult_c <-  RunUMAP(srat_Adult_c,dims=seq(30)) 
srat_Adult_c <-  FindNeighbors(srat_Adult_c,dims=seq(30)) 
srat_Adult_c <-  FindClusters(srat_Adult_c,resolution=1) 
Adult_c_DR <-  as.data.frame(srat_Adult_c@reductions$umap@cell.embeddings)
colnames(Adult_c_DR) = c('RD1','RD2')
Adult_c_DR$Cluster = factor(srat_Adult_c@meta.data[rownames(Adult_c_DR),'RNA_snn_res.1'])

#Adult_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
Adult_c_DR$sftpc = Adult_c.sc$toc["Sftpc", rownames(Adult_c_DR)]
Adult_c_DR$scgb1a1 = Adult_c.sc$toc["Scgb1a1", rownames(Adult_c_DR)]
Adult_c_DR$ptprc = Adult_c.sc$toc["Ptprc", rownames(Adult_c_DR)]
Adult_c_DR$pdgfra = Adult_c.sc$toc["Pdgfra", rownames(Adult_c_DR)]
Adult_c_DR$pecam1 = Adult_c.sc$toc["Pecam1", rownames(Adult_c_DR)]
Adult_c_DR$dcn = Adult_c.sc$toc["Dcn", rownames(Adult_c_DR)]
ggplot(Adult_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(Adult_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(Adult_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(Adult_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(Adult_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(Adult_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(Adult_c.sc, "Scgb1a1", Adult_c_DR)
plotMarkerMap(Adult_c.sc, "Sftpc", Adult_c_DR)
plotMarkerMap(Adult_c.sc, "Ptprc", Adult_c_DR)
plotMarkerMap(Adult_c.sc, "Pdgfra", Adult_c_DR)
plotMarkerMap(Adult_c.sc, "Pecam1", Adult_c_DR)
plotMarkerMap(Adult_c.sc, "Dcn", Adult_c_DR)

#Adult_c, which genes are expressed most highly in the background? 
head(Adult_c.sc$soupProfile[order(Adult_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(Adult_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_Adult_c = estimateNonExpressingCells(Adult_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(Adult_c_DR$Cluster, rownames(Adult_c_DR)))
plotMarkerMap(Adult_c.sc, geneSet = background_RNA_genes, DR = Adult_c_DR, useToEst = useToEst_Adult_c)

## calculate contamination and adjust the counts. 
Adult_c.sc = setClusters(Adult_c.sc, Adult_c_DR$Cluster)
Adult_c.sc = calculateContaminationFraction(Adult_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_Adult_c)
head(Adult_c.sc$metaData)

## adjust the counts based on contamination fraction
out_Adult_c = adjustCounts(Adult_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(Adult_c.sc$toc > 0)
cntStrained = rowSums(out_Adult_c > 0)
mostZeroed_Adult_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_Adult_c
plotChangeMap(Adult_c.sc, out_Adult_c, "Cd74", Adult_c_DR)
plotChangeMap(Adult_c.sc, out_Adult_c, "Scgb3a2", Adult_c_DR)
plotChangeMap(Adult_c.sc, out_Adult_c, "Scgb1a1", Adult_c_DR)
plotChangeMap(Adult_c.sc, out_Adult_c, "Sftpc", Adult_c_DR)
plotChangeMap(Adult_c.sc, out_Adult_c, "Dcn", Adult_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./Adult_c_counts", out_Adult_c)

###############################################################

#Adult_2_c non-epithelial specific sort SoupX Pipeline

#Read matrix files into R for Adult_2_c
Adult_2_c.raw.matrix <- Read10X(data.dir = "./raw_feature_bc_matrix JS-11 P64/", gene.column = 2, unique.features = TRUE)
Adult_2_c.filtered.matrix <- Read10X(data.dir = "./filtered_feature_bc_matrix JS11 P64-2", gene.column = 2, unique.features = TRUE)

#Data come out as a list of two matrices (i.e. 'Gene Expression' and 'Antibody Capture')
#Subset the list of matrices to get the gene expression matrix 
Adult_2_c.raw.matrix.GE <- Adult_2_c.raw.matrix$`Gene Expression`
Adult_2_c.filtered.matrix.GE <- Adult_2_c.filtered.matrix$`Gene Expression`

#Manually create SoupChannel
Adult_2_c.sc <- SoupChannel(tod = Adult_2_c.raw.matrix.GE, toc = Adult_2_c.filtered.matrix.GE, metaData = NULL, soupRange = c(0, 10), keepDroplets = FALSE)

# tSNE coordinates for Adult_2_c data
srat_Adult_2_c <-  CreateSeuratObject(Adult_2_c.sc$toc)
srat_Adult_2_c <-  NormalizeData(srat_Adult_2_c)
srat_Adult_2_c <-  FindVariableFeatures(srat_Adult_2_c) 
srat_Adult_2_c <-  ScaleData(srat_Adult_2_c) 
srat_Adult_2_c <-  RunPCA(srat_Adult_2_c,pcs.compute=30) 
srat_Adult_2_c <-  RunUMAP(srat_Adult_2_c,dims=seq(30)) 
srat_Adult_2_c <-  FindNeighbors(srat_Adult_2_c,dims=seq(30)) 
srat_Adult_2_c <-  FindClusters(srat_Adult_2_c,resolution=1) 
Adult_2_c_DR <-  as.data.frame(srat_Adult_2_c@reductions$umap@cell.embeddings)
colnames(Adult_2_c_DR) = c('RD1','RD2')
Adult_2_c_DR$Cluster = factor(srat_Adult_2_c@meta.data[rownames(Adult_2_c_DR),'RNA_snn_res.1'])

#Adult_2_c Plotting of marker maps to check for RNA contamination. See if the gene of interest is expessed in different clusters
Adult_2_c_DR$sftpc = Adult_2_c.sc$toc["Sftpc", rownames(Adult_2_c_DR)]
Adult_2_c_DR$scgb1a1 = Adult_2_c.sc$toc["Scgb1a1", rownames(Adult_2_c_DR)]
Adult_2_c_DR$ptprc = Adult_2_c.sc$toc["Ptprc", rownames(Adult_2_c_DR)]
Adult_2_c_DR$pdgfra = Adult_2_c.sc$toc["Pdgfra", rownames(Adult_2_c_DR)]
Adult_2_c_DR$pecam1 = Adult_2_c.sc$toc["Pecam1", rownames(Adult_2_c_DR)]
Adult_2_c_DR$dcn = Adult_2_c.sc$toc["Dcn", rownames(Adult_2_c_DR)]
ggplot(Adult_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = sftpc > 0))
ggplot(Adult_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = scgb1a1 > 0))
ggplot(Adult_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = ptprc > 0))
ggplot(Adult_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pdgfra > 0))
ggplot(Adult_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = pecam1 > 0))
ggplot(Adult_2_c_DR, aes(RD1, RD2)) + geom_point(aes(colour = dcn > 0))
plotMarkerMap(Adult_2_c.sc, "Scgb1a1", Adult_2_c_DR)
plotMarkerMap(Adult_2_c.sc, "Sftpc", Adult_2_c_DR)
plotMarkerMap(Adult_2_c.sc, "Ptprc", Adult_2_c_DR)
plotMarkerMap(Adult_2_c.sc, "Pdgfra", Adult_2_c_DR)
plotMarkerMap(Adult_2_c.sc, "Pecam1", Adult_2_c_DR)
plotMarkerMap(Adult_2_c.sc, "Dcn", Adult_2_c_DR)

#Adult_2_c, which genes are expressed most highly in the background? 
head(Adult_2_c.sc$soupProfile[order(Adult_2_c.sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution plots the distribution of the observed to expected expression for marker genes
plotMarkerDistribution(Adult_2_c.sc)

#Specify background RNA genes for non-epithelial cell specific sort
background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1")

#Specify background RNA genes for epithelial specific cell sort
#background_RNA_genes = c("Dcn", "Bgn", "Aspn", "Ecm2", "Fos", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Lyz1", "Lyz2", "Mgp", "Postn", "Scgb1a1", "Sftpc", "Hopx", "Ager", "Krt19", "Cldn4", "Foxj1", "Krt5", "Sfn", "Pecam1")

useToEst_Adult_2_c = estimateNonExpressingCells(Adult_2_c.sc, nonExpressedGeneList = list(SLRP = background_RNA_genes), clusters = setNames(Adult_2_c_DR$Cluster, rownames(Adult_2_c_DR)))
plotMarkerMap(Adult_2_c.sc, geneSet = background_RNA_genes, DR = Adult_2_c_DR, useToEst = useToEst_Adult_2_c)

## calculate contamination and adjust the counts. 
Adult_2_c.sc = setClusters(Adult_2_c.sc, Adult_2_c_DR$Cluster)
Adult_2_c.sc = calculateContaminationFraction(Adult_2_c.sc, list(SLRP = background_RNA_genes), useToEst = useToEst_Adult_2_c)
head(Adult_2_c.sc$metaData)

## adjust the counts based on contamination fraction
out_Adult_2_c = adjustCounts(Adult_2_c.sc)

## rechecking some plots and genes to check the genes are mostly specific to  clusters. 
cntSoggy = rowSums(Adult_2_c.sc$toc > 0)
cntStrained = rowSums(out_Adult_2_c > 0)
mostZeroed_Adult_2_c = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed_Adult_2_c
plotChangeMap(Adult_2_c.sc, out_Adult_2_c, "Cd74", Adult_2_c_DR)
plotChangeMap(Adult_2_c.sc, out_Adult_2_c, "Scgb3a2", Adult_2_c_DR)
plotChangeMap(Adult_2_c.sc, out_Adult_2_c, "Scgb1a1", Adult_2_c_DR)
plotChangeMap(Adult_2_c.sc, out_Adult_2_c, "Sftpc", Adult_2_c_DR)
plotChangeMap(Adult_2_c.sc, out_Adult_2_c, "Dcn", Adult_2_c_DR)

### write the file to a new count matrix. saves similar to original one.
DropletUtils:::write10xCounts("./Adult_2_c_counts", out_Adult_2_c)

###############################################################
