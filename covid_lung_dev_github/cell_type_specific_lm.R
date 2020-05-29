#load in needed libraries
library(Seurat)
library(dplyr)
#Set working direcotry (probably folder where you are keeping the data for the project)
setwd("~/Desktop/covid_python/")
#Set random seed to get reproducible results
set.seed(33)

#N.B. the location where the anovas are generated is hard coded
anova_generator_epi <- function(seuratObject, gene) {
  location = toString("~/Desktop/covid_python/anova_generator")
  if(!file.exists(location)) {
    dir.create(location)
  }
  
  #Get gene expression data from Seurat object
  Stats_table <- data.frame(FetchData(object = seuratObject, vars = c('celltype', 'timepoint', str_replace(gene, "\\.", "-"))))
  #Pull celltypes from Stats table
  Stats_table_AT1 <- Stats_table[Stats_table$celltype == "AT1",]
  Stats_table_AT2 <- Stats_table[Stats_table$celltype == "AT2",]
  Stats_table_Secretory <- Stats_table[Stats_table$celltype == "Secretory",]
  Stats_table_Ciliated <- Stats_table[Stats_table$celltype == "Ciliated",]
  Stats_table_PNEC <- Stats_table[Stats_table$celltype == "PNEC",]
  
  ct <- c("AT1", "AT2", "Secretory", "Ciliated", "PNEC")
  
  gene_anova <- data.frame()
  for (i in 1:length(ct)) {
    lm_i <- lm(eval(as.symbol(gene)) ~ timepoint, data = eval(as.symbol(paste("Stats_table_", ct[i], sep = ""))))
    anova_i <- anova(lm_i)
    anova_i$celltype <- ct[i]
    gene_anova <- rbind(gene_anova, anova_i)
  }
  write.csv(gene_anova, toString(paste(location, "/", gene, "_anova.csv", sep = "")))
}

COVID_SO <- readRDS("epi_lung_dev_covid.rds")
Idents(COVID_SO) <- "leiden"
COVID_SO$celltype <- Idents(COVID_SO)

anova_generator_epi(COVID_SO, "Ace2")
anova_generator_epi(COVID_SO, "Tmprss2")
anova_generator_epi(COVID_SO, "Furin")
anova_generator_epi(COVID_SO, "Ctsb")
anova_generator_epi(COVID_SO, "Bsg")
anova_generator_epi(COVID_SO, "Irf3")
anova_generator_epi(COVID_SO, "Tnf")
anova_generator_epi(COVID_SO, "Mapk1")
anova_generator_epi(COVID_SO, "Mapk3")
anova_generator_epi(COVID_SO, "Nfkb1")
anova_generator_epi(COVID_SO, "Nfkb2")
anova_generator_epi(COVID_SO, "Cxcl15")
anova_generator_epi(COVID_SO, "Ifnar1")
anova_generator_epi(COVID_SO, "Ifnar2")
anova_generator_epi(COVID_SO, "Ifitm1")
anova_generator_epi(COVID_SO, "Ifitm3")
anova_generator_epi(COVID_SO, "Jak1")
anova_generator_epi(COVID_SO, "Mx1")
anova_generator_epi(COVID_SO, "H2.K1")
anova_generator_epi(COVID_SO, "H2.Aa")
anova_generator_epi(COVID_SO, "H2.Ab1")
anova_generator_epi(COVID_SO, "Cd74")
anova_generator_epi(COVID_SO, "H2.DMb1")