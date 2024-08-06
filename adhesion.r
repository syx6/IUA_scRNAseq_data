library(Seurat)
library(patchwork)
library(DoubletFinder)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(DoubletFinder)
library(Seurat)
library(tidyverse)
library("BuenColors")
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(AnnotationHub) 
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(tidyverse)
library(reshape2)
library(circlize)
library(ggsignif)

fileNames <- list.files(pattern="adhesion_") 
n = length(fileNames)

for (s in 1:n){
  sample = fileNames[s]
  print(sample)
  assign(sample, Read10X(data.dir = sample))
  
}

fileNames <- list.files(pattern="control") 
n = length(fileNames)

for (s in 1:n){
  sample = fileNames[s]
  print(sample)
  assign(sample, Read10X(data.dir = sample))
  
}


adhesion_1 <- CreateSeuratObject(counts = adhesion_1, project = "adhesion1", min.cells = 3, min.features = 300)
adhesion_2 <- CreateSeuratObject(counts = adhesion_2, project = "adhesion2", min.cells = 3, min.features = 300)
adhesion_3 <- CreateSeuratObject(counts = adhesion_3, project = "adhesion3", min.cells = 3, min.features = 300)
adhesion_4 <- CreateSeuratObject(counts = adhesion_4, project = "adhesion4", min.cells = 3, min.features = 300)

control1 <- CreateSeuratObject(counts = control1, project = "control1", min.cells = 3, min.features = 300)
control2 <- CreateSeuratObject(counts = control2, project = "control2", min.cells = 3, min.features = 300)
control3 <- CreateSeuratObject(counts = control3, project = "control3", min.cells = 3, min.features = 300)
control4 <- CreateSeuratObject(counts = control4, project = "control4", min.cells = 3, min.features = 300)
control5 <- CreateSeuratObject(counts = control3, project = "control5", min.cells = 3, min.features = 300)
control6 <- CreateSeuratObject(counts = control4, project = "control6", min.cells = 3, min.features = 300)
control7 <- CreateSeuratObject(counts = control3, project = "control7", min.cells = 3, min.features = 300)



adhesion_1[["percent.mt"]] <- PercentageFeatureSet(adhesion_1, pattern = "^MT-")
adhesion_2[["percent.mt"]] <- PercentageFeatureSet(adhesion_2, pattern = "^MT-")
adhesion_3[["percent.mt"]] <- PercentageFeatureSet(adhesion_3, pattern = "^MT-")
adhesion_4[["percent.mt"]] <- PercentageFeatureSet(adhesion_4, pattern = "^MT-")
control1[["percent.mt"]] <- PercentageFeatureSet(control1, pattern = "^MT-")
control2[["percent.mt"]] <- PercentageFeatureSet(control2, pattern = "^MT-")
control3[["percent.mt"]] <- PercentageFeatureSet(control3, pattern = "^MT-")
control4[["percent.mt"]] <- PercentageFeatureSet(control4, pattern = "^MT-")


VlnPlot(adhesion_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)
VlnPlot(adhesion_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)
VlnPlot(adhesion_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)
VlnPlot(adhesion_4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)
VlnPlot(control1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)
VlnPlot(control2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)
VlnPlot(control3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)
VlnPlot(control4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.25)

adhesion_1 <- subset(adhesion_1, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)
adhesion_2 <- subset(adhesion_2, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)
adhesion_3 <- subset(adhesion_3, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)
adhesion_4 <- subset(adhesion_4, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)

control1 <- subset(control1, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)
control2 <- subset(control2, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)
control3 <- subset(control3, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)
control4 <- subset(control4, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 20)


control1 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> control1

ElbowPlot(control1, ndims=50)


control2 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> control2

ElbowPlot(control2, ndims=50)

#################
control3 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> control3

ElbowPlot(control3, ndims=50)


control4 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> control4

ElbowPlot(control4, ndims=50)

#####################

sweep.res.list_SM <- paramSweep_v3(control1, PCs = 1:40,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- control1@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(control1)*8*1e-6
nExp_poi <- round(DoubletRate*length(control1@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
control1.DF <- doubletFinder_v3(control1, PCs = 1:40, pN = pN_value, 
                                  pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


sweep.res.list_SM <- paramSweep_v3(control2, PCs = 1:40,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- control2@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(control2)*8*1e-6
nExp_poi <- round(DoubletRate*length(control2@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
control2.DF <- doubletFinder_v3(control2, PCs = 1:40, pN = pN_value, 
                                pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


sweep.res.list_SM <- paramSweep_v3(control3, PCs = 1:30,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- control3@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(control3)*8*1e-6
nExp_poi <- round(DoubletRate*length(control3@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
control3.DF <- doubletFinder_v3(control3, PCs = 1:30, pN = pN_value, 
                                pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

sweep.res.list_SM <- paramSweep_v3(control4, PCs = 1:30,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- control4@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(control4)*8*1e-6
nExp_poi <- round(DoubletRate*length(control4@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
control4.DF <- doubletFinder_v3(control4, PCs = 1:30, pN = pN_value, 
                                pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)



adhesion_1 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> adhesion_1

ElbowPlot(adhesion_1, ndims=50)

adhesion_2 %>% Seurat::RunUMAP( dims = 1:30, verbose = FALSE) %>% 
  Seurat::RunTSNE( dims = 1:30, verbose = FALSE) %>%
  Seurat::FindNeighbors( dims = 1:30, verbose = FALSE)%>% 
  Seurat::FindClusters(verbose = FALSE) -> adhesion_2

adhesion_2 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> adhesion_2

ElbowPlot(adhesion_2, ndims=50)

adhesion_2 %>% Seurat::RunUMAP( dims = 1:30, verbose = FALSE) %>% 
  Seurat::RunTSNE( dims = 1:30, verbose = FALSE) %>%
  Seurat::FindNeighbors( dims = 1:30, verbose = FALSE)%>% 
  Seurat::FindClusters(verbose = FALSE) -> adhesion_2

adhesion_3 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> adhesion_3

ElbowPlot(adhesion_3, ndims=50)

adhesion_3 %>% Seurat::RunUMAP( dims = 1:30, verbose = FALSE) %>% 
  Seurat::RunTSNE( dims = 1:30, verbose = FALSE) %>%
  Seurat::FindNeighbors( dims = 1:30, verbose = FALSE)%>% 
  Seurat::FindClusters(verbose = FALSE) -> adhesion_3

adhesion_4 %>%  NormalizeData()%>% 
  FindVariableFeatures()%>% 
  ScaleData()  %>% 
  Seurat::RunPCA( verbose = FALSE) -> adhesion_4

ElbowPlot(adhesion_4, ndims=50)

adhesion_4 %>% Seurat::RunUMAP( dims = 1:30, verbose = FALSE) %>% 
  Seurat::RunTSNE( dims = 1:30, verbose = FALSE) %>%
  Seurat::FindNeighbors( dims = 1:30, verbose = FALSE)%>% 
  Seurat::FindClusters(verbose = FALSE) -> adhesion_4

sweep.res.list_SM <- paramSweep_v3(adhesion_4, PCs = 1:30,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- adhesion_4@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(adhesion_4)*8*1e-6
nExp_poi <- round(DoubletRate*length(adhesion_4@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
adhesion_4.DF <- doubletFinder_v3(adhesion_4, PCs = 1:30, pN = pN_value, 
                             pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


sweep.res.list_SM <- paramSweep_v3(adhesion_3, PCs = 1:30,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- adhesion_3@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(adhesion_3)*8*1e-6
nExp_poi <- round(DoubletRate*length(adhesion_3@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
adhesion_3.DF <- doubletFinder_v3(adhesion_3, PCs = 1:30, pN = pN_value, 
                                  pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

sweep.res.list_SM <- paramSweep_v3(adhesion_2, PCs = 1:30,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- adhesion_2@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(adhesion_2)*8*1e-6
nExp_poi <- round(DoubletRate*length(adhesion_2@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
adhesion_2.DF <- doubletFinder_v3(adhesion_2, PCs = 1:30, pN = pN_value, 
                                  pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


sweep.res.list_SM <- paramSweep_v3(adhesion_1, PCs = 1:30,sct = FALSE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- adhesion_1@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)   
DoubletRate = ncol(adhesion_1)*8*1e-6
nExp_poi <- round(DoubletRate*length(adhesion_1@meta.data$orig.ident))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
adhesion_1.DF <- doubletFinder_v3(adhesion_1, PCs = 1:30, pN = pN_value, 
                                  pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


adhesion_1.DF@meta.data$Doublet <- adhesion_1.DF@meta.data$DF.classifications_0.25_0.29_145
adhesion_1.DF.Singlet <- subset(x=adhesion_1.DF, subset = Doublet =='Singlet')

adhesion_2.DF@meta.data$Doublet <- adhesion_2.DF@meta.data$DF.classifications_0.25_0.3_323
adhesion_2.DF.Singlet <- subset(x=adhesion_2.DF, subset = Doublet =='Singlet')

adhesion_3.DF@meta.data$Doublet <- adhesion_3.DF@meta.data$DF.classifications_0.25_0.29_128
adhesion_3.DF.Singlet <- subset(x=adhesion_3.DF, subset = Doublet =='Singlet')

adhesion_4.DF@meta.data$Doublet <- adhesion_4.DF@meta.data$DF.classifications_0.25_0.27_172
adhesion_4.DF.Singlet <- subset(x=adhesion_4.DF, subset = Doublet =='Singlet')

control1.DF@meta.data$Doublet <- control1.DF@meta.data$DF.classifications_0.25_0.3_590
control1.DF.Singlet <- subset(x=control1.DF, subset = Doublet =='Singlet')

control2.DF@meta.data$Doublet <- control2.DF@meta.data$DF.classifications_0.25_0.22_941
control2.DF.Singlet <- subset(x=control2.DF, subset = Doublet =='Singlet')

control3.DF@meta.data$Doublet <- control3.DF@meta.data$DF.classifications_0.25_0.3_109
control3.DF.Singlet <- subset(x=control3.DF, subset = Doublet =='Singlet')

control4.DF@meta.data$Doublet <- control4.DF@meta.data$DF.classifications_0.25_0.05_111
control4.DF.Singlet <- subset(x=control4.DF, subset = Doublet =='Singlet')



DimPlot(adhesion_4.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()
DimPlot(adhesion_3.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()
DimPlot(adhesion_2.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()
DimPlot(adhesion_1.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()

DimPlot(control2.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()
DimPlot(control1.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()

DimPlot(control3.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()
DimPlot(control4.DF.Singlet, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()



adhesion.total.Singlet <- merge(adhesion_1.DF.Singlet, y = c(adhesion_2.DF.Singlet,adhesion_3.DF.Singlet,
                                                             adhesion_4.DF.Singlet,control1.DF.Singlet,
                                                             control2.DF.Singlet),
                            add.cell.ids = c("adhesion_1", "adhesion_2",'adhesion_3',
                                             'adhesion_4','control1','control2'), project = "adhesion")


adhesion.total.Singlet.new <- merge(control3.DF.Singlet, y = c(control4.DF.Singlet,adhesion.total.Singlet),
                                    add.cell.ids = c('control3','control4',''), project = "adhesion")




adhesion.total.Singlet.SCT <- SCTransform(adhesion.total.Singlet, 
                                           vars.to.regress = "percent.mt", 
                                           verbose = FALSE) 

adhesion.total.Singlet.SCT <- RunPCA(adhesion.total.Singlet.SCT, verbose = FALSE)
ElbowPlot(adhesion.total.Singlet.SCT, ndims=50)
adhesion.total.Singlet.SCT <- RunUMAP(adhesion.total.Singlet.SCT, dims = 1:50, verbose = FALSE)
adhesion.total.Singlet.SCT <- RunTSNE(adhesion.total.Singlet.SCT, dims = 1:50, verbose = FALSE)
adhesion.total.Singlet.SCT <- FindNeighbors(adhesion.total.Singlet.SCT, dims = 1:50)
adhesion.total.Singlet.SCT <- FindClusters(adhesion.total.Singlet.SCT)

Idents(adhesion.total.Singlet.SCT) <- 'SCT_snn_res.0.8'

DimPlot(adhesion.total.Singlet.SCT, reduction = "umap",label=TRUE) + ggsci::scale_color_igv() + theme_test()

FeaturePlot(adhesion.total.Singlet.SCT, features = c('KRT8','COL1A1',
                                                     'PTPRC','ACTA2'))


adhesion.total.Singlet.SCT@meta.data <- adhesion.total.Singlet.SCT@meta.data %>%
  mutate(celltype=case_when(      SCT_snn_res.0.8 == "6" |
                                  SCT_snn_res.0.8 == "18" |
                                  SCT_snn_res.0.8 == "25" ~ 'epi',
                                  SCT_snn_res.0.8 == "1" |
                                  SCT_snn_res.0.8 == "2" |
                                  SCT_snn_res.0.8 == "3" |
                                  SCT_snn_res.0.8 == "9" |
                                  SCT_snn_res.0.8 == "12" |
                                  SCT_snn_res.0.8 == "16" |
                                  SCT_snn_res.0.8 == "15" |
                                  SCT_snn_res.0.8 == "17" |
                                  SCT_snn_res.0.8 == "29" |
                                  SCT_snn_res.0.8 == "21" |
                                  SCT_snn_res.0.8 == "31" ~ 'fib',
                                  SCT_snn_res.0.8 == "0" |
                                  SCT_snn_res.0.8 == "5" |
                                  SCT_snn_res.0.8 == "4" |
                                  SCT_snn_res.0.8 == "11" |
                                  SCT_snn_res.0.8 == "13" |
                                  SCT_snn_res.0.8 == "20" |
                                  SCT_snn_res.0.8 == "26" |
                                  SCT_snn_res.0.8 == "27" |
                                  SCT_snn_res.0.8 == "10" |
                                  SCT_snn_res.0.8 == "7" |
                                  SCT_snn_res.0.8 == "8" |
                                  SCT_snn_res.0.8 == "14" |
                                  SCT_snn_res.0.8 == "19" |
                                  SCT_snn_res.0.8 == "30" |
                                  SCT_snn_res.0.8 == "23" | 
                                  SCT_snn_res.0.8 == "28" |
                                  SCT_snn_res.0.8 == "24"  |
                                  SCT_snn_res.0.8 == "22" | 
                                  SCT_snn_res.0.8 == "32"  ~ 'immune cells',
                                TRUE ~ as.character(SCT_snn_res.0.8)))


adhesion.total.Singlet.SCT@meta.data <- adhesion.total.Singlet.SCT@meta.data %>%
  mutate(group=case_when(           orig.ident == "adhesion1" |
                                    orig.ident == "adhesion2" |
                                    orig.ident == "adhesion3" |
                                    orig.ident == "adhesion4" ~ 'adhesion',
                                    orig.ident == "control1" |
                                    orig.ident == "control2" ~ 'normal',
                                  TRUE ~ as.character(orig.ident)))


adhesion.total.list <- SplitObject(adhesion.total.Singlet, split.by = "orig.ident")

adhesion.total.list <- lapply(X = adhesion.total.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adhesion.total.list)
adhesion.total.list <- PrepSCTIntegration(object.list = adhesion.total.list, anchor.features = features)
adhesion.total.anchors <- FindIntegrationAnchors(object.list = adhesion.total.list, normalization.method = "SCT", 
                                             anchor.features = features)
adhesion.total.sct <- IntegrateData(anchorset = adhesion.total.anchors, normalization.method = "SCT")


adhesion.total.sct <- RunPCA(adhesion.total.sct, verbose = FALSE)

ElbowPlot(adhesion.total.sct, ndims=50)

adhesion.total.sct <- RunUMAP(adhesion.total.sct, dims = 1:40, verbose = FALSE)
adhesion.total.sct <- RunTSNE(adhesion.total.sct, dims = 1:40, verbose = FALSE)
adhesion.total.sct <- FindNeighbors(adhesion.total.sct, dims = 1:40)
adhesion.total.sct <- FindClusters(adhesion.total.sct,resolution = 0.8)

DimPlot(adhesion.total.sct, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()

saveRDS(adhesion.total.sct, file='adhesion.total.sct.rds')
saveRDS(adhesion.total.Singlet,file='adhesion.total.Singlet.rds')

saveRDS(adhesion.total.Singlet.SCT, file='adhesion.total.Singlet.SCT.rds')
saveRDS(adhesion.total.Singlet,file='adhesion.total.Singlet.rds')

unciliated <- c('KRT18')
Macrophage <- c('CD68', 'MS4A4A', 'MS4A7')
Fibro <- c('COL1A1', 'COL3A1', 'COL1A2','DCN','COL6A1')
Endo <- c('AQP1', 'MYCT1', 'CDH5', 'PECAM1')
T <- c('CD3D','CD8A')
Mast <- c('TPSB2', 'TPSAB1')
NK <- c('KLRF1','NCAM1','GNLY','NKG7')

B <- c('CD19','CD79B','MS4A1','CD79A')
NKT <- c('CD3D','CD3G','GZMB')
Pla <- c('MZB1','JCHAIN')
Myo <- c('ACTA2','MYH11')
Ciliated <- c('SNTN','FOXJ1','CDHR3')
Fib <- c('LUM')
MONO <- c('S100A8')
DC <- c('CD1C')


DotPlot(adhesion.total.sct, features = c(T,unciliated,Macrophage,NK,Ciliated,Fib,MONO,Mast,B,Endo,DC,Pla), 
        cols =  c("white", "deeppink3"),scale.by='size') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_text(size=12,family="Arial")) +
  theme(axis.text.x=element_text(size=12,angle=90,family="Arial"),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.title.x=element_text(size = 18,family="Arial"),
        axis.title.y=element_text(size = 18,family="Arial")) + ylab('Celltype') + xlab('Marker gene') +
  theme(legend.key.height = unit(10,'pt'))


DimPlot(adhesion.total.Singlet.SCT, reduction = "umap",label = TRUE) + ggsci::scale_color_igv() + theme_test()

FeaturePlot(adhesion.total.Singlet.SCT, features = c('CD68'))

adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.8 == "16" ~ 'Endothelium',
                            integrated_snn_res.0.8 == "9" ~ 'Myofibroblast',
                            integrated_snn_res.0.8 == "19" ~ 'Mast cells',
                            integrated_snn_res.0.8 == "10" |
                              integrated_snn_res.0.8 == "12"|
                              integrated_snn_res.0.8 == "21"|
                              integrated_snn_res.0.8 == "20"~ 'Macrophage',
                            integrated_snn_res.0.8 == "0"|
                              integrated_snn_res.0.8 == "3"|
                              integrated_snn_res.0.8 == "5"~ 'Stromal fibroblast',
                            integrated_snn_res.0.8 == "18"  ~ 'B',
                            integrated_snn_res.0.8 == "24"  ~ 'Pla',
                            integrated_snn_res.0.8 == "14"  ~ 'Ciliated',
                            integrated_snn_res.0.8 == "6"|
                              integrated_snn_res.0.8 == "13"|
                              integrated_snn_res.0.8 == "11"~ 'Unciliated',
                            integrated_snn_res.0.8 == "2"|
                              integrated_snn_res.0.8 == "7"|
                              integrated_snn_res.0.8 == "15"|
                              integrated_snn_res.0.8 == "17"~ 'T',
                            integrated_snn_res.0.8 == "1"|
                              integrated_snn_res.0.8 == "4"|
                              integrated_snn_res.0.8 == "23"~ 'NK',
                            TRUE ~ as.character(integrated_snn_res.0.8)
  ))


cluster_markers_all <- adhesion.total.sct@misc$cluster_markers_all <- Seurat::FindAllMarkers(object = adhesion.total.sct, 
                                                                                                assay = "RNA",
                                                                                                slot = "data",
                                                                                                verbose = TRUE, 
                                                                                                only.pos = TRUE, 
                                                                                                logfc.threshold = 0.9,
                                                                                                min.pct = 0.5)


adhesion.total.sct.new <- subset(adhesion.total.sct,idents = c('T','NK','Macrophage',
                                                               'Unciliated','Stromal fibroblast','Mast cells','Ciliated',
                                                               'Endothelium','Myofibroblast','Pla','B'))


saveRDS(adhesion.total.sct.new, file='adhesion.total.sct.new.rds')


#绘图

DimPlot(adhesion.total.sct.new,reduction = 'umap') + 
  ggsci::scale_color_igv() + guides(color=guide_legend(override.aes = list(size=6))) + theme_classic() + 
  theme(legend.text=element_text(size=12))


adhesion.total.DF.Singlet.average <- AverageExpression( adhesion.total.DF.Singlet, 
                                                    return.seurat=TRUE )


DoHeatmap(adhesion.total.DF.Singlet.average, features = unique(cluster_markers_all$gene), angle = 45,
          size = 4, draw.lines = FALSE) + 
  scale_fill_gradientn(colors = c('#104E8B', 'white', '#EE1289')) +  
  theme(legend.title=element_text(size=14,family="serif")) +
  theme(legend.text = element_text(colour = 'black', size = 12,,family="serif")) +
  theme(legend.key.size=unit(0.5,'cm')) + theme(legend.key.width=unit(0.5,'cm')) 


VlnPlot(adhesion.total.DF.Singlet, features = c("nFeature_RNA"), ncol = 1, pt.size=0,cols=colorRampPalette(brewer.pal(8,'Accent'))(10)) +
  labs(title="Detected gene per cell type") + theme(title=element_text(size=16,face="bold")) + xlab(NULL)



DimPlot(adhesion.total.sct.new,reduction = 'umap',label=TRUE,pt.size = 0.3) + 
  scale_colour_manual(values =  c('#3300FF','#FFCCFF','#FF9966','#FF6600','#CCFF00','#CCCCFF',
                                  '#CCCCCC','#66FF00','#00EE00','#EE0000','#66CCFF'))+
  guides(color=guide_legend(override.aes = list(size=6))) 


adhesion.total.sct.new@meta.data <- adhesion.total.sct.new@meta.data %>% mutate(group=case_when(orig.ident == "adhesion1" |
                                                                      orig.ident == "adhesion2" |
                                                                      orig.ident == "adhesion3" |
                                                                      orig.ident == "adhesion4" ~ "IU",
                                                                    orig.ident == "control1" |
                                                                      orig.ident == "control2"  ~ "Normal"))


# 细胞 注释GLI+

FibToMyo@meta.data <- FibToMyo@meta.data %>% mutate(CellID=rownames(FibToMyo@meta.data))
GLI1_annote <- GLI1@meta.data %>% dplyr::select(CellType, CellID)
FibToMyo@meta.data <- left_join(FibToMyo@meta.data, GLI1_annote, by = "CellID")
row.names(FibToMyo@meta.data) <- FibToMyo@meta.data$CellID

FibToMyo@meta.data <- FibToMyo@meta.data %>% mutate(CellType=case_when(CellType == 'GL1+' & celltype == 'Myofibroblast' ~ 'GL1+ Myofibroblast',
                                                                       CellType == 'GL1+' & celltype == 'Stromal fibroblast' ~ 'GL1+ Stromal fibroblast',
                                                                       TRUE ~ celltype))

FibToMyo@meta.data <- FibToMyo@meta.data %>% mutate(state=case_when(CellType == 'Myofibroblast' ~ 'GL1-',
                                                                    CellType == 'Stromal fibroblast' ~ 'GL1-',
                                                                    TRUE ~ 'GL1+'))


### fibroblast-to-myofibroblast

FibToMyo <- subset(adhesion.total.sct.new,idents = c('Stromal fibroblast','Myofibroblast'))
Fib <- subset(FibToMyo,idents = c('0','1','3'))

Idents(FibToMyo) <- 'integrated_snn_res.0.8'

DefaultAssay(FibToMyo) <- "integrated"
FibToMyo <- RunPCA(FibToMyo, verbose = FALSE)
ElbowPlot(FibToMyo, ndims=50)
FibToMyo <- RunUMAP(FibToMyo, dims = 1:30, verbose = FALSE)
FibToMyo <- RunTSNE(FibToMyo, dims = 1:30, verbose = FALSE)
FibToMyo <- FindNeighbors(FibToMyo, dims = 1:30)
FibToMyo <- FindClusters(FibToMyo, resolution = 0.2)

#Fib <- subset(FibToMyo,idents = c('Stromal fibroblast','GL1+ Stromal fibroblast'))
#Epi <- subset(adhesion.total.sct.new,idents = c('Ciliated','Unciliated'))


#FibToMyo@meta.data <- FibToMyo@meta.data %>% mutate(group=case_when(orig.ident == "adhesion1" |
 #                                     orig.ident == "adhesion2" |
  #                                    orig.ident == "adhesion3" |
   #                                   orig.ident == "adhesion4" ~ "IU",
    #                                orig.ident == "control1" |
     #                                 orig.ident == "control2"  ~ "Normal"))


#subset(x = FibToMyo, subset = IHH > 0, slot = 'counts')

#GLI <- subset(FibToMyo,idents = c('GL1+ Stromal fibroblast', 'GL1+ Myofibroblast'))
#Epi <- subset(adhesion.total.sct.new,idents = c('Ciliated','Unciliated'))
#Epi_IHH <- subset(x = Epi, subset = IHH > 0, slot = 'counts')
EMT <- subset(adhesion.total.sct.new,idents = c('Ciliated','Unciliated','Stromal fibroblast','Myofibroblast'))



##拟时分析

data <- as(as.matrix(FibToMyo@assays$RNA@counts), 'sparseMatrix')
pd <- FibToMyo@meta.data
pdd <- pd %>% mutate(Type=case_when(orig.ident == "adhesion1" |
                                      orig.ident == "adhesion2" |
                                      orig.ident == "adhesion3" |
                                      orig.ident == "adhesion4" ~ "IU",
                                      orig.ident == "control1" |
                                      orig.ident == "control2"  ~ "Normal"))
pddd <- new('AnnotatedDataFrame', data = pdd)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,phenoData = pddd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],fullModelFormulaStr = "~integrated_snn_res.0.2")
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1800]

#ordering_genes <- row.names (subset(diff_test_res, qval < 0.000000001))

monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)

plot_cell_trajectory(monocle_cds, color_by = "celltype",show_branch_points = F,cell_size = 2) +
  scale_color_manual(values= c("#0073C2FF","pink"))  + guides(color=guide_legend(override.aes = list(size=8)))

plot_cell_trajectory(adhesion_monocle_cds, color_by = "celltype",show_branch_points = F,cell_size = 0.9)+ 
  guides(color=guide_legend(override.aes = list(size=4)))+ facet_wrap(~Type, ncol = 2) 


###
cluster_markers_all_FibToMyo <- FibToMyo@misc$cluster_markers_all <- Seurat::FindAllMarkers(object = FibToMyo, 
                                                                                            assay = "RNA",
                                                                                             slot = "data",
                                                                                             verbose = TRUE, 
                                                                                             only.pos = TRUE, 
                                                                                             logfc.threshold = 0.5,
                                                                                             min.pct = 0.5)


cluster_markers_all_Fib <- Fib@misc$cluster_markers_all <- Seurat::FindAllMarkers(object = Fib, 
                                                                                            assay = "RNA",
                                                                                            slot = "data",
                                                                                            verbose = TRUE, 
                                                                                            only.pos = TRUE, 
                                                                                            logfc.threshold = 0.5,
                                                                                            min.pct = 0.5)


###
Idents(adhesion.total.sct.new) <- 'celltype'
adhesion.total.sct.new.fibro <- subset(adhesion.total.sct.new,idents = "Stromal fibroblast")
Idents(adhesion.total.sct.new.fibro) <- 'group'
fibro_diff_all <- adhesion.total.sct.new.fibro@misc$fibro_diff_all <- 
  FindMarkers(adhesion.total.sct.new.fibro,
              ident.1 = c('IU'), 
              ident.2 = c('Normal'), 
              verbose = FALSE,
              logfc.threshold = 0)


############# 整合public data

adhesion.total.Singlet <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/endometirum_adhesion/adhesion.total.Singlet.rds')



adhesion.total.list <- SplitObject(adhesion.total.Singlet.new, split.by = "orig.ident")

adhesion.total.list <- lapply(X = adhesion.total.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adhesion.total.list)
adhesion.total.list <- PrepSCTIntegration(object.list = adhesion.total.list, anchor.features = features)
adhesion.total.anchors <- FindIntegrationAnchors(object.list = adhesion.total.list, normalization.method = "SCT", 
                                                 anchor.features = features)
adhesion.total.sct <- IntegrateData(anchorset = adhesion.total.anchors, normalization.method = "SCT")


adhesion.total.sct <- RunPCA(adhesion.total.sct, verbose = FALSE)

ElbowPlot(adhesion.total.sct, ndims=50)

adhesion.total.sct <- RunUMAP(adhesion.total.sct, dims = 1:40, verbose = FALSE)
adhesion.total.sct <- RunTSNE(adhesion.total.sct, dims = 1:40, verbose = FALSE)
adhesion.total.sct <- FindNeighbors(adhesion.total.sct, dims = 1:40)
adhesion.total.sct <- FindClusters(adhesion.total.sct,resolution = 1.2)

DimPlot(adhesion.total.sct, reduction = "umap",label=T) + ggsci::scale_color_igv() + theme_test()



#######
saveRDS(adhesion.total.Singlet.new, 
        file='/szrmyy/wangjgLab/scRNA/xiasy/endometirum_adhesion/adhesion.total.Singlet.new.rds')

saveRDS(adhesion.total.sct, 
        file='/szrmyy/wangjgLab/scRNA/xiasy/endometirum_adhesion/adhesion.total.sct.rds')

adhesion.total.sct <- readRDS(file='/szrmyy/wangjgLab/scRNA/xiasy/endometirum_adhesion/adhesion.total.sct.rds')


##########细胞注释

DotPlot(adhesion.total.sct, features = c(T,unciliated,Macrophage,NK,Ciliated,Fibro,Mast,Endo,Pla,B), 
        cols =  c("white", "deeppink3"),scale.by='size') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_text(size=12,family="serif",face ='bold')) +
  theme(axis.text.x=element_text(size=12,angle=90,family="serif", face ='bold'),
        axis.text.y=element_text(size=12,family="serif",face ='bold'),
        axis.title.x=element_text(size = 18,family="serif",face ='bold'),
        axis.title.y=element_text(size = 18,family="serif",face ='bold')) + ylab('Celltype') + xlab('Marker gene') +
  theme(legend.key.height = unit(10,'pt'))


FeaturePlot(adhesion.total.sct, features = c('COL1A1','ACTA2'),label=TRUE)


adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.8 == "19" ~ 'Endothelium',
                            integrated_snn_res.0.8 == "1" |
                              integrated_snn_res.0.8 == "3" |
                              integrated_snn_res.0.8 == "5" |
                              integrated_snn_res.0.8 == "10" | 
                              integrated_snn_res.0.8 == "15" | 
                              integrated_snn_res.0.8 == "25" ~ 'Fibroblast',
                            integrated_snn_res.0.8 == "17"  ~ 'Ciliated',
                            integrated_snn_res.0.8 == "4"|
                              integrated_snn_res.0.8 == "7"|
                              integrated_snn_res.0.8 == "9"|
                              integrated_snn_res.0.8 == "11"|
                              integrated_snn_res.0.8 == "18"~ 'Unciliated',
                            integrated_snn_res.0.8 == "0"|
                              integrated_snn_res.0.8 == "6"  ~ 'NK',
                            integrated_snn_res.0.8 == "2"|
                              integrated_snn_res.0.8 == "8"|
                              integrated_snn_res.0.8 == "21"|
                              integrated_snn_res.0.8 == "22"|
                              integrated_snn_res.0.8 == "24" ~ 'T',
                            integrated_snn_res.0.8 == "20"~ 'Mast',
                            TRUE ~ as.character(integrated_snn_res.1.2)
  ))


adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>%
  mutate(CellType=case_when(
    celltype == "34"~ 'Pla',
    celltype == "25"~ 'B',
    celltype == "16" |
    celltype == "19" ~ 'Macrophage',
    celltype == "20"~ 'NK',
    celltype == "33"~ 'Dendritic',
    celltype == "5" |
      celltype == "7" |
      celltype == "13" |
      celltype == "14" |
      celltype == "15"  ~ 'T',
    celltype == "18"  ~ 'Fibroblast',
    celltype == "10" |
      celltype == "12"  ~ 'Unciliated',
    celltype == "9"  ~ ' Neutrophil',
    
    TRUE ~ as.character(celltype)
  ))


adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>%
  mutate(group=case_when(           orig.ident == "adhesion1" |
                                      orig.ident == "adhesion2" |
                                      orig.ident == "adhesion3" |
                                      orig.ident == "adhesion4" ~ 'adhesion',
                                    orig.ident == "control1" |
                                      orig.ident == "control2" |
                                      orig.ident == "control3" |
                                      orig.ident == "control4"~ 'normal',
                                    TRUE ~ as.character(orig.ident)))


cluster_markers_all <- adhesion.total.sct@misc$cluster_markers_all <- Seurat::FindAllMarkers(object = adhesion.total.sct, 
                                                                                             assay = "RNA",
                                                                                             slot = "data",
                                                                                             verbose = TRUE, 
                                                                                             only.pos = TRUE, 
                                                                                             logfc.threshold = 0.9,
                                                                                             min.pct = 0.5)


adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>% dplyr::select("orig.ident","nCount_RNA","nFeature_RNA",
                                                                               "percent.mt","Doublet", "nCount_SCT",
                                                                               "nFeature_SCT","integrated_snn_res.0.8",
                                                                               "CellType","group")

adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>%
  mutate(CellType=case_when(CellType == "T" ~ 'T',
                            CellType == "Unciliated" ~ 'UNCILIA',
                            CellType == "Macrophage" ~ 'MAC',
                            CellType == "NK" ~ 'NK',
                            CellType == "Ciliated" ~ 'CILIA',
                            CellType == "Fibroblast" ~ 'FIB',
                            CellType == " Neutrophil" ~ 'MONO',
                            CellType == "Mast" ~ 'MAST',
                            CellType == "B" ~ 'B',
                            CellType == "Endothelium" ~ 'ENDO',
                            CellType == "Dendritic" ~ 'DC',
                            CellType == "Pla" ~ 'Pla',
                            TRUE ~ as.character(CellType)))

### 绘图 Figure1

DimPlot(adhesion.total.sct,reduction = 'tsne',label=TRUE,pt.size = 0.3) + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 

######### "CD3D"   "CD8A"   "KRT18"  "CD68"   "MS4A4A" "GNLY"   "NKG7"   "SNTN"   
##### "FOXJ1"  "LUM"    "S100A8" "TPSB2"  "TPSAB1" "MS4A1"  "CD79A"  "CDH5"  
######## "PECAM1" "CD1C"   "MZB1"   "JCHAIN"

DotPlot(adhesion.total.sct, features = c(T,unciliated,Macrophage,NK,Ciliated,Fib,MONO,Mast,B,Endo,DC,Pla), 
        cols =  c("white", "deeppink3"),scale.by='size') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_text(size=12,family="Arial")) +
  theme(axis.text.x=element_text(size=12,angle=90,family="Arial"),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.title.x=element_text(size = 18,family="Arial"),
        axis.title.y=element_text(size = 18,family="Arial")) + ylab('Celltype') + xlab('Marker gene') +
  theme(legend.key.height = unit(10,'pt'))

FeaturePlot(adhesion.total.sct, features = c('PTPRC'),reduction = 'tsne',label=TRUE)

adhesion.total.sct.average <- AverageExpression(adhesion.total.sct, return.seurat=TRUE )


DoHeatmap(adhesion.total.sct.average, features = unique(top10$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#00EE00','#FFCCFF','#FF9966','#CCCCFF',
                                                  '#66CCFF','#EE0000')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink1')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 

top1 <- adhesion.total.sct@misc$cluster_markers_all %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

#######统计细胞数
cell_num <- data.frame(table(Idents(adhesion.total.sct)))

ggplot(data=cell_num,mapping=aes(x=Var1,y=Freq,fill=Var1))+ geom_bar(stat="identity") + coord_flip() + 
  scale_fill_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  theme_classic()+
  theme(axis.text.x=element_text(size=11,family="Arial"),
        axis.text.y=element_text(size=11,family="Arial"),
        axis.title.x=element_text(size = 11,family="Arial"),
        axis.title.y=element_text(size = 12,family="Arial")) + ylab('Cell number')

###########细胞比例 整体
library(cowplot)
library(ggplot2)
library(dplyr)
library(ggpubr)

Idents(adhesion.total.sct) <- 'group'
adhesion <- subset(adhesion.total.sct, idents='adhesion')
normal <- subset(adhesion.total.sct, idents='normal')

Idents(adhesion) <- 'CellType'
Idents(normal) <- 'CellType'

cellprop <- rbind(as.data.frame(prop.table(table(Idents(adhesion), adhesion$orig.ident))),
                  as.data.frame(prop.table(table(Idents(normal), normal$orig.ident))))

colnames(cellprop)<-c("celltype","group","proportion")

cellprop$group <- as.character(cellprop$group)

cellprop <- cellprop %>% 
  mutate(group1=case_when(group == "adhesion1" ~ 'adhesion',
                          group == "adhesion2"  ~ 'adhesion',
                          group == "adhesion3" ~ 'adhesion',
                          group == "adhesion4" ~ 'adhesion',
                          group == "control1" ~ 'normal',
                          group == "control2" ~ 'normal',
                          group == "control3" ~ 'normal',
                          group == "control4" ~ 'normal',
                          TRUE ~ '1'))

ggplot(cellprop %>% filter_all(any_vars(str_detect(., pattern="Fibroblast"))), aes(celltype, weight = proportion, fill = group1)) +
  geom_bar(color = "black", width = .7, position = 'dodge') + theme_classic2() + coord_flip() + scale_fill_manual(values = c('pink','lightblue'))

library(reshape2)
cellprop <- dcast(cellprop,group~celltype, value.var = "proportion")
rownames(cellprop) <- cellprop[,1]
cellprop <- cellprop[,-1]

###添加分组信息
sample <- c("adhesion1","adhesion2","adhesion3","adhesion4","control1","control2","control3","control4")
group <- c("adhesion","adhesion","adhesion","adhesion","normal","normal","normal",'normal')
samples <- data.frame(sample, group)#创建数据框

rownames(samples)=samples$sample
cellprop$sample <- samples[rownames(cellprop),'sample']#R添加列
cellprop$group <- samples[rownames(cellprop),'group']

## cellprop.fibro <- cellprop %>% dplyr::select('SFRP4+ Myofibroblast', 'SCARA5+ Fibroblast', 'RGS5+ Myofibroblast', 'CCN1+ Fibroblast',sample,group)

pplist = list()
sce_groups = colnames(cellprop)[1:12]


for(group_ in sce_groups){
  cellper_  = cellprop %>% dplyr::select(one_of(c('sample','group',group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25,size=2) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10),legend.position = 'none') + 
    labs(title = group_,y='Percentage') + scale_fill_manual(values = c('deeppink','lightblue'))+
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  pplist[[group_]] = pp1
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("normal", "adhesion") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
  
}

plot_grid(pplist[['T']],
          pplist[['NK']],
          pplist[['MAC']],
          pplist[['DC']],
          pplist[['UNCILIA']],
          pplist[['FIB']],
          pplist[['MAST']],
          pplist[['ENDO']],
          pplist[['B']],
          pplist[['Pla']],
          pplist[['CILIA']],
          pplist[['MONO']]
          )
#####################
library(Seurat)
library(ggplot2)
# 构建函数
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
                          ...) {
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


StackedVlnPlot(adhesion.total.sct, top1$gene, pt.size=0, cols= c('lightblue','#FFCCFF','#FF9966',
                                                                                           'pink','#CCFF00','#CCCCFF',
                                                                                           '#CCCCCC','mediumpurple','#66FF00',
                                                                                           '#EE0000','#66CCFF','tomato')) + coord_flip()





####差异基因

Idents(adhesion.total.sct) <- 'CellType'
adhesion.fibro <- subset(adhesion.total.sct,idents = "FIB")
Idents(adhesion.fibro) <- 'group'
adhesion.fibro_diff_all <- adhesion.fibro@misc$adhesion.fibro_diff_all <- 
  FindMarkers(adhesion.fibro,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.T <- subset(adhesion.total.sct,idents = "T")
Idents(adhesion.T) <- 'group'
adhesion.T_diff_all <- adhesion.T@misc$adhesion.T_diff_all <- 
  FindMarkers(adhesion.T,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Unciliated <- subset(adhesion.total.sct,idents = "Unciliated")
Idents(adhesion.Unciliated) <- 'group'
adhesion.Unciliated_diff_all <- adhesion.Unciliated@misc$adhesion.Unciliated_diff_all <- 
  FindMarkers(adhesion.Unciliated,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Macrophage <- subset(adhesion.total.sct,idents = "Macrophage")
Idents(adhesion.Macrophage) <- 'group'
adhesion.Macrophage_diff_all <- adhesion.Macrophage@misc$adhesion.Macrophage_diff_all <- 
  FindMarkers(adhesion.Macrophage,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.NK <- subset(adhesion.total.sct,idents = "NK")
Idents(adhesion.NK) <- 'group'
adhesion.NK_diff_all <- adhesion.NK@misc$adhesion.NK_diff_all <- 
  FindMarkers(adhesion.NK,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Ciliated <- subset(adhesion.total.sct,idents = "Ciliated")
Idents(adhesion.Ciliated) <- 'group'
adhesion.Ciliated_diff_all <- adhesion.Ciliated@misc$adhesion.Ciliated_diff_all <- 
  FindMarkers(adhesion.Ciliated,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Neutrophil <- subset(adhesion.total.sct,idents = " Neutrophil")
Idents(adhesion.Neutrophil) <- 'group'
adhesion.Neutrophil_diff_all <- adhesion.Neutrophil@misc$adhesion.Neutrophil_diff_all <- 
  FindMarkers(adhesion.Neutrophil,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Mast <- subset(adhesion.total.sct,idents = "Mast")
Idents(adhesion.Mast) <- 'group'
adhesion.Mast_diff_all <- adhesion.Mast@misc$adhesion.Mast_diff_all <- 
  FindMarkers(adhesion.Mast,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.B <- subset(adhesion.total.sct,idents = "B")
Idents(adhesion.B) <- 'group'
adhesion.B_diff_all <- adhesion.B@misc$adhesion.B_diff_all <- 
  FindMarkers(adhesion.B,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Endothelium <- subset(adhesion.total.sct,idents = "Endothelium")
Idents(adhesion.Endothelium) <- 'group'
adhesion.Endothelium_diff_all <- adhesion.Endothelium@misc$adhesion.Endothelium_diff_all <- 
  FindMarkers(adhesion.Endothelium,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Dendritic <- subset(adhesion.total.sct,idents = "Dendritic")
Idents(adhesion.Dendritic) <- 'group'
adhesion.Dendritic_diff_all <- adhesion.Dendritic@misc$adhesion.Dendritic_diff_all <- 
  FindMarkers(adhesion.Dendritic,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion.total.sct) <- 'CellType'
adhesion.Pla <- subset(adhesion.total.sct,idents = "Pla")
Idents(adhesion.Pla) <- 'group'
adhesion.Pla_diff_all <- adhesion.Pla@misc$adhesion.Pla_diff_all <- 
  FindMarkers(adhesion.Pla,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)


###差异基因list


adhesion.fibro_diff_all <- adhesion.fibro_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Fibroblast') %>% 
  mutate(gene=row.names(adhesion.fibro_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.T_diff_all <- adhesion.T_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='T') %>% 
  mutate(gene=row.names(adhesion.T_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Unciliated_diff_all <- adhesion.Unciliated_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Unciliated') %>% 
  mutate(gene=row.names(adhesion.Unciliated_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Macrophage_diff_all <- adhesion.Macrophage_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Macrophage') %>% 
  mutate(gene=row.names(adhesion.Macrophage_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.NK_diff_all <- adhesion.NK_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='NK') %>% 
  mutate(gene=row.names(adhesion.NK_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Ciliated_diff_all <- adhesion.Ciliated_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Ciliated') %>% 
  mutate(gene=row.names(adhesion.Ciliated_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Neutrophil_diff_all <- adhesion.Neutrophil_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Neutrophil') %>% 
  mutate(gene=row.names(adhesion.Neutrophil_diff_all)) %>% 
  rownames_to_column("gene1")


adhesion.Mast_diff_all <- adhesion.Mast_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Mast') %>% 
  mutate(gene=row.names(adhesion.Mast_diff_all)) %>% 
  rownames_to_column("gene1")


adhesion.B_diff_all <- adhesion.B_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='B') %>% 
  mutate(gene=row.names(adhesion.B_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Endothelium_diff_all <- adhesion.Endothelium_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Endothelium') %>% 
  mutate(gene=row.names(adhesion.Endothelium_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Dendritic_diff_all <- adhesion.Dendritic_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Dendritic') %>% 
  mutate(gene=row.names(adhesion.Dendritic_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Pla_diff_all <- adhesion.Pla_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Pla') %>% 
  mutate(gene=row.names(adhesion.Pla_diff_all)) %>% 
  rownames_to_column("gene1")


adhesion.diff.gene <- rbind(adhesion.B_diff_all,adhesion.Ciliated_diff_all,adhesion.Dendritic_diff_all,
                            adhesion.Endothelium_diff_all,adhesion.Macrophage_diff_all,adhesion.Mast_diff_all,
                            adhesion.Neutrophil_diff_all, adhesion.NK_diff_all, adhesion.Pla_diff_all,
                            adhesion.T_diff_all, adhesion.Unciliated_diff_all, adhesion.fibro_diff_all)


saveRDS(adhesion.diff.gene, file='adhesion.diff.gene.rds')
saveRDS(adhesion.total.sct, file='adhesion.total.sct.rds')


adhesion.diff.gene <- adhesion.diff.gene %>%
  mutate(celltype=case_when(celltype == "T" ~ 'T',
                            celltype == "Unciliated" ~ 'UNCILIA',
                            celltype == "Macrophage" ~ 'MAC',
                            celltype == "NK" ~ 'NK',
                            celltype == "Ciliated" ~ 'CILIA',
                            celltype == "Fibroblast" ~ 'FIB',
                            celltype == " Neutrophil" ~ 'MONO',
                            celltype == "Mast" ~ 'MAST',
                            celltype == "B" ~ 'B',
                            celltype == "Endothelium" ~ 'ENDO',
                            celltype == "Dendritic" ~ 'DC',
                            celltype == "Pla" ~ 'Pla',
                            TRUE ~ as.character(celltype)))


################################Epi

adhesion_epi <- subset(adhesion.total.sct,idents = c("UNCILIA","CILIA"))

##########################################
adhesion_epi.list <- SplitObject(adhesion_epi, split.by = "orig.ident")

adhesion_epi.list <- lapply(X = adhesion_epi.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adhesion_epi.list)
adhesion_epi.list <- PrepSCTIntegration(object.list = adhesion_epi.list, anchor.features = features)
adhesion_epi.anchors <- FindIntegrationAnchors(object.list = adhesion_epi.list, normalization.method = "SCT", 
                                                 anchor.features = features)
adhesion_epi.sct <- IntegrateData(anchorset = adhesion_epi.anchors, normalization.method = "SCT")

#adhesion_epi <- SCTransform(adhesion_epi, vars.to.regress = "percent.mt",verbose = FALSE)

DefaultAssay(adhesion_epi) <- 'integrated'
adhesion_epi <- RunPCA(adhesion_epi, verbose = FALSE)
ElbowPlot(adhesion_epi, ndims=50)
adhesion_epi <- RunUMAP(adhesion_epi, dims = 1:30, verbose = FALSE)
adhesion_epi <- RunTSNE(adhesion_epi, dims = 1:30, verbose = FALSE)
adhesion_epi <- FindNeighbors(adhesion_epi, dims = 1:30)
adhesion_epi <- FindClusters(adhesion_epi,resolution = 0.8)

Idents(adhesion_epi) <- 'CellType'

DefaultAssay(adhesion_epi) <- 'RNA'
FeaturePlot(adhesion_epi, features = c('EPCAM','UCA1','SNTN','FOXJ1','KLF5','WFDC2','UCA1','COL1A1'),reduction = 'umap',label=TRUE)

DotPlot(adhesion_epi, features = c('EPCAM','UCA1','SNTN','FOXJ1','KLF5','WFDC2','TACSTD2'), 
        cols =  c("white", "deeppink3"),scale.by='size') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_text(size=12,family="Arial")) +
  theme(axis.text.x=element_text(size=12,angle=90,family="Arial"),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.title.x=element_text(size = 18,family="Arial"),
        axis.title.y=element_text(size = 18,family="Arial")) + ylab('Celltype') + xlab('Marker gene') +
  theme(legend.key.height = unit(10,'pt'))




Idents(adhesion_epi.sct) <- 'integrated_snn_res.0.8'
DimPlot(adhesion_epi.sct,reduction = 'umap',label=TRUE,pt.size = 0.3,split.by = 'group')+
  guides(color=guide_legend(override.aes = list(size=6))) 

DimPlot(adhesion_epi,reduction = 'umap',label=TRUE,pt.size = 0.3) + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 


adhesion_epi@meta.data <- adhesion_epi@meta.data %>%
  mutate(CellType=case_when(SCT_snn_res.0.8 == "2" |
                            SCT_snn_res.0.8 == "17" |
                            SCT_snn_res.0.8 == "18"  ~ 'UNCILIA',
                            SCT_snn_res.0.8 == "4" |
                            SCT_snn_res.0.8 == "7" |
                            SCT_snn_res.0.8 == "9" |
                            SCT_snn_res.0.8 == "10" |
                            SCT_snn_res.0.8 == "11" |
                            SCT_snn_res.0.8 == "21" |
                            SCT_snn_res.0.8 == "19"  ~ 'FIB',
                            TRUE ~ as.character('CILIA')
  ))

adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>% mutate(CellID=rownames(adhesion.total.sct@meta.data))
adhesion_epi@meta.data <- adhesion_epi@meta.data %>% mutate(CellID=rownames(adhesion_epi@meta.data))

adhesion_epi_annote <- adhesion_epi@meta.data %>% dplyr::select(CellType, CellID)
adhesion.total.sct@meta.data <- left_join(adhesion.total.sct@meta.data, 
                                          adhesion_epi_annote, by = "CellID")

adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>%
  mutate(CellType=case_when(CellType.y == 'CILIA' ~ 'UNCILIA',
                            CellType.y == 'UNCILIA' ~ 'CILIA',
                            CellType.y == 'FIB' ~ 'FIB',
                            TRUE ~ as.character(CellType.x)))
  
adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>% dplyr::select(-CellType.y,-CellType.x)
rownames(adhesion.total.sct@meta.data) <- adhesion.total.sct@meta.data$CellID
  
  
########上皮细胞注释：
DefaultAssay(adhesion_epi) <- 'integrated'
adhesion_epi <- RunPCA(adhesion_epi, verbose = FALSE)
ElbowPlot(adhesion_epi, ndims=50)
adhesion_epi <- RunUMAP(adhesion_epi, dims = 1:30, verbose = FALSE)
adhesion_epi <- RunTSNE(adhesion_epi, dims = 1:30, verbose = FALSE)
adhesion_epi <- FindNeighbors(adhesion_epi, dims = 1:30)
adhesion_epi <- FindClusters(adhesion_epi,resolution = 2)


adhesion_epi@meta.data <- adhesion_epi@meta.data %>%
  mutate(CellType=case_when(SCT_snn_res.2 == "2" |
                            SCT_snn_res.2 == "13" |
                            SCT_snn_res.2 == "19" |
                            SCT_snn_res.2 == "20" |
                            SCT_snn_res.2 == "25" |
                            SCT_snn_res.2 == "29" ~ 'CILIA',
                            SCT_snn_res.2 == "0" |
                              SCT_snn_res.2 == "1" |
                              SCT_snn_res.2 == "4" |
                              SCT_snn_res.2 == "5" |
                              SCT_snn_res.2 == "6" |
                              SCT_snn_res.2 == "7" |
                              SCT_snn_res.2 == "9" |
                              SCT_snn_res.2 == "24" |
                              SCT_snn_res.2 == "15" |
                              SCT_snn_res.2 == "22" |
                              SCT_snn_res.2 == "27" ~ 'G_UNCILIA',
                            TRUE ~ as.character('L_UNCILIA')))
  


#############Stromal cells
adhesion_stromal <- subset(adhesion.total.sct,idents = c("UNCILIA","CILIA",
                                                         'FIB','ENDO'))


########stromal细胞注释：

adhesion_stromal.list <- SplitObject(adhesion_stromal, split.by = "orig.ident")

adhesion_stromal.list <- lapply(X = adhesion_stromal.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adhesion_stromal.list)
adhesion_stromal.list <- PrepSCTIntegration(object.list = adhesion_stromal.list, anchor.features = features)
adhesion_stromal.anchors <- FindIntegrationAnchors(object.list = adhesion_stromal.list, normalization.method = "SCT", 
                                               anchor.features = features)
adhesion_stromal.sct <- IntegrateData(anchorset = adhesion_stromal.anchors, normalization.method = "SCT")




DefaultAssay(adhesion_stromal) <- 'integrated'
adhesion_stromal.sct <- RunPCA(adhesion_stromal.sct, verbose = FALSE)
ElbowPlot(adhesion_stromal, ndims=50)
adhesion_stromal.sct <- RunUMAP(adhesion_stromal.sct, dims = 1:40, verbose = FALSE)
adhesion_stromal.sct <- RunTSNE(adhesion_stromal.sct, dims = 1:40, verbose = FALSE)
adhesion_stromal.sct <- FindNeighbors(adhesion_stromal.sct, dims = 1:40)
adhesion_stromal.sct <- FindClusters(adhesion_stromal.sct,resolution = 2)

DefaultAssay(adhesion_stromal.sct) <- 'RNA'
Idents(adhesion_stromal) <- 'CellType'
Idents(adhesion_stromal) <- 'integrated_snn_res.2'

FeaturePlot(adhesion_stromal, features = 'MCAM')
FeaturePlot(adhesion_stromal.sct, features = c('MCAM','CLDN5','KRT18','COL1A1','EPCAM'),label=T)

saveRDS(adhesion_stromal.sct ,file='~/data/adhesion/adhesion_stromal.sct.rds')

DimPlot(adhesion_stromal.sct,reduction = 'umap',label=TRUE,pt.size = 0.3) + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 


Idents(adhesion_stromal.sct) <- 'integrated_snn_res.2'
Idents(adhesion_stromal.sct) <- 'CellType'
DotPlot(adhesion_stromal.sct, features = c('MCAM'), 
        cols =  c("white", "deeppink3"),scale.by='size') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_text(size=12,family="Arial")) +
  theme(axis.text.x=element_text(size=12,angle=90,family="Arial"),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.title.x=element_text(size = 18,family="Arial"),
        axis.title.y=element_text(size = 18,family="Arial")) + ylab('Celltype') + xlab('Marker gene') +
  theme(legend.key.height = unit(10,'pt'))

FeaturePlot(adhesion_stromal.sct, features = c('UCA1','EPCAM','COL1A1'))

adhesion_stromal.sct@meta.data <- adhesion_stromal.sct@meta.data %>%
  mutate(CellType=case_when(integrated_snn_res.2 == "26" |
                            integrated_snn_res.2 == "7" |
                            integrated_snn_res.2 == "30" ~ 'SMC',
                            TRUE ~ as.character(CellType)))
  

########diff gene fib smo

Idents(adhesion_stromal.sct) <- 'CellType'
adhesion.Fib <- subset(adhesion_stromal.sct,idents = "FIB")
Idents(adhesion.Fib) <- 'group'
adhesion.Fib_diff_all <- adhesion.Fib@misc$adhesion.Fib_diff_all <- 
  FindMarkers(adhesion.Fib,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

Idents(adhesion_stromal.sct) <- 'CellType'
adhesion.smc <- subset(adhesion_stromal.sct,idents = "SMC")
Idents(adhesion.smc) <- 'group'
adhesion.smc_diff_all <- adhesion.smc@misc$adhesion.smc_diff_all <- 
  FindMarkers(adhesion.smc,
              ident.1 = c('adhesion'), 
              ident.2 = c('normal'), 
              verbose = FALSE,
              logfc.threshold = 0)

adhesion.smc_diff_all <- adhesion.smc_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='SMC') %>% 
  mutate(gene=row.names(adhesion.smc_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.Fib_diff_all <- adhesion.Fib_diff_all %>% 
  mutate(sig=case_when(avg_logFC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_logFC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_logFC > -0.25 | avg_logFC < 0.25 ~ 'non')) %>%
  mutate(celltype='Fibroblast') %>% 
  mutate(gene=row.names(adhesion.Fib_diff_all)) %>% 
  rownames_to_column("gene1")

adhesion.diff.gene <- adhesion.diff.gene %>% filter(celltype != 'Fibroblast')
adhesion.diff.gene <- rbind(adhesion.diff.gene, adhesion.smc_diff_all, adhesion.Fib_diff_all)
saveRDS(adhesion.diff.gene, file='adhesion.diff.gene.rds')

##############Fibroblast

#############Fibroblast
adhesion_Fibroblast <- subset(adhesion_stromal.sct,idents = c('FIB'))

########stromal细胞注释：

DefaultAssay(adhesion_Fibroblast.sct) <- 'integrated'
adhesion_Fibroblast.sct <- RunPCA(adhesion_Fibroblast.sct, verbose = FALSE)
ElbowPlot(adhesion_Fibroblast.sct, ndims=50)
adhesion_Fibroblast.sct <- RunUMAP(adhesion_Fibroblast.sct, dims = 1:50, verbose = FALSE)
adhesion_Fibroblast.sct <- RunTSNE(adhesion_Fibroblast.sct, dims = 1:50, verbose = FALSE)
adhesion_Fibroblast.sct <- FindNeighbors(adhesion_Fibroblast.sct, dims = 1:50)
adhesion_Fibroblast.sct <- FindClusters(adhesion_Fibroblast.sct,resolution = 2)
##############Macrophage

adhesion_Macrophage <- subset(adhesion.total.sct,idents = c('MAC','MONO'))

########Macrophage细胞注释：

adhesion_Macrophage.list <- SplitObject(adhesion_Macrophage, split.by = "orig.ident")

adhesion_Macrophage.list <- lapply(X = adhesion_Macrophage.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adhesion_Macrophage.list)
adhesion_Macrophage.list <- PrepSCTIntegration(object.list = adhesion_Macrophage.list, anchor.features = features)
adhesion_Macrophage.anchors <- FindIntegrationAnchors(object.list = adhesion_Macrophage.list, normalization.method = "SCT", 
                                                   anchor.features = features)
adhesion_Macrophage.sct <- IntegrateData(anchorset = adhesion_Macrophage.anchors, normalization.method = "SCT")



adhesion_Macrophage.sct <- RunPCA(adhesion_Macrophage.sct, verbose = FALSE)
ElbowPlot(adhesion_Macrophage.sct, ndims=50)
adhesion_Macrophage.sct <- RunUMAP(adhesion_Macrophage.sct, dims = 1:40, verbose = FALSE)
adhesion_Macrophage.sct <- RunTSNE(adhesion_Macrophage.sct, dims = 1:40, verbose = FALSE)
adhesion_Macrophage.sct <- FindNeighbors(adhesion_Macrophage.sct, dims = 1:40)
adhesion_Macrophage.sct <- FindClusters(adhesion_Macrophage.sct,resolution = 2)

DimPlot(adhesion_Macrophage.sct,reduction = 'umap',label=TRUE,pt.size = 0.3) + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 

##############T cells
adhesion_T <- subset(adhesion.total.sct,idents = c('T'))

########Macrophage细胞注释：

adhesion_T.list <- SplitObject(adhesion_T, split.by = "orig.ident")

adhesion_T.list <- lapply(X = adhesion_T.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adhesion_T.list)
adhesion_T.list <- PrepSCTIntegration(object.list = adhesion_T.list, anchor.features = features)
adhesion_T.list.anchors <- FindIntegrationAnchors(object.list = adhesion_T.list, normalization.method = "SCT", 
                                                      anchor.features = features)
adhesion_T.list.sct <- IntegrateData(anchorset = adhesion_T.list.anchors, normalization.method = "SCT")


adhesion_T.sct <- RunPCA(adhesion_T.list.sct, verbose = FALSE)
ElbowPlot(adhesion_T.sct, ndims=50)
adhesion_T.sct <- RunUMAP(adhesion_T.sct, dims = 1:40, verbose = FALSE)
adhesion_T.sct <- RunTSNE(adhesion_T.sct, dims = 1:40, verbose = FALSE)
adhesion_T.sct <- FindNeighbors(adhesion_T.sct, dims = 1:40)
adhesion_T.sct <- FindClusters(adhesion_T.sct,resolution = 2)


DimPlot(adhesion_T.sct,reduction = 'umap',label=TRUE,pt.size = 0.3) + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 



saveRDS(adhesion.diff.gene, file='adhesion.diff.gene.rds')
saveRDS(adhesion_stromal.sct, file='adhesion_stromal.sct.rds')
saveRDS(adhesion_Macrophage.sct, file='adhesion_Macrophage.sct.rds')
saveRDS(adhesion_Fibroblast.sct, file='adhesion_Fibroblast.sct.rds')
saveRDS(adhesion_T.sct, file='adhesion_T.sct.rds')


############### marker gene epi
##   SOX9, LGR5: SOX9+ basal
##   FOXJ,PIFO: Ciliated
##   MUC12, CDC20B,CCNO, HES6: Pre-ciliated
###  PTGS1: Lumenal
###  SCGB2A2: Glandular


adhesion.total.sct <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion.total.sct.rds')
adhesion.diff.gene <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion.diff.gene.rds')


######### total go terms

library(scRNAtoolVis)


adhesion.diff.gene <- adhesion.diff.gene %>% mutate(avg_log2FC = avg_logFC) %>% 
  mutate(cluster = celltype)

adheson.diff.immunegene <- adhesion.diff.gene %>% dplyr::filter(celltype == 'B'|
                                                                  celltype == 'Dendritic'|
                                                                  celltype == 'Macrophage'|
                                                                  celltype == 'Mast'|
                                                                  celltype == 'Neutrophil'|
                                                                  celltype == 'NK'|
                                                                  celltype == 'Pla'|
                                                                  celltype == 'T'
                                                                  )


jjVolcano(diffData = adheson.diff.immunegene,
          legend.position=c(0.95,0.05),
          log2FC.cutoff = 0.25,
          col.type = "updown",
          adjustP.cutoff=0.05,
          topGeneN=3,
          size  = 3.5,
          fontface = 'italic',
          flip = T)


Fibroblast.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Fibroblast')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

SMC.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'SMC')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Unciliated.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Unciliated')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

T.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'T')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Pla.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Pla')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


NK.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'NK')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Neutrophil.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Neutrophil')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Mast.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Mast')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Macrophage.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Macrophage')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Endothelium.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Endothelium')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Dendritic.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Dendritic')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Ciliated.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Ciliated')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

B.gene.up.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'B')  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


####### down go terms

Fibroblast.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Fibroblast')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

SMC.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'SMC')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Unciliated.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Unciliated')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

T.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'T')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Pla.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Pla')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


NK.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'NK')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Neutrophil.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Neutrophil')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Mast.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Mast')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Macrophage.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Macrophage')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Endothelium.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Endothelium')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Dendritic.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Dendritic')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

Ciliated.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'Ciliated')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

B.gene.down.ego <- adhesion.diff.gene %>%
  dplyr::filter(celltype == 'B')  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


###################################################
### go terms merge
########## downgene
B.gene.down.ego.result <- B.gene.down.ego@result %>% mutate(celltype='B') %>% mutate(sig='down')
Ciliated.gene.down.ego.result <- Ciliated.gene.down.ego@result %>% mutate(celltype='Ciliated') %>% mutate(sig='down')
Dendritic.gene.down.ego.result <- Dendritic.gene.down.ego@result %>% mutate(celltype='Dendritic') %>% mutate(sig='down')
Endothelium.gene.down.ego.result <- Endothelium.gene.down.ego@result %>% mutate(celltype='Endothelium') %>% mutate(sig='down')
Macrophage.gene.down.ego.result <- Macrophage.gene.down.ego@result %>% mutate(celltype='Macrophage') %>% mutate(sig='down')
Mast.gene.down.ego.result <- Mast.gene.down.ego@result %>% mutate(celltype='Mast') %>% mutate(sig='down')
Neutrophil.gene.down.ego.result <- Neutrophil.gene.down.ego@result %>% mutate(celltype='Neutrophil') %>% mutate(sig='down')
NK.gene.down.ego.result <- NK.gene.down.ego@result %>% mutate(celltype='NK') %>% mutate(sig='down')
Pla.gene.down.ego.result <- Pla.gene.down.ego@result %>% mutate(celltype='Pla') %>% mutate(sig='down')
T.gene.down.ego.result <- T.gene.down.ego@result %>% mutate(celltype='T') %>% mutate(sig='down')
Unciliated.gene.down.ego.result <- Unciliated.gene.down.ego@result %>% mutate(celltype='Unciliated') %>% mutate(sig='down')
SMC.gene.down.ego.result <- SMC.gene.down.ego@result %>% mutate(celltype='SMC') %>% mutate(sig='down')
Fibroblast.gene.down.ego.result <- Fibroblast.gene.down.ego@result %>% mutate(celltype='Fibroblast') %>% mutate(sig='down')

### go terms merge
########## upgene
B.gene.up.ego.result <- B.gene.up.ego@result %>% mutate(celltype='B') %>% mutate(sig='up')
Ciliated.gene.up.ego.result <- Ciliated.gene.up.ego@result %>% mutate(celltype='Ciliated') %>% mutate(sig='up')
Dendritic.gene.up.ego.result <- Dendritic.gene.up.ego@result %>% mutate(celltype='Dendritic') %>% mutate(sig='up')
Endothelium.gene.up.ego.result <- Endothelium.gene.up.ego@result %>% mutate(celltype='Endothelium') %>% mutate(sig='up')
Macrophage.gene.up.ego.result <- Macrophage.gene.up.ego@result %>% mutate(celltype='Macrophage') %>% mutate(sig='up')
Mast.gene.up.ego.result <- Mast.gene.up.ego@result %>% mutate(celltype='Mast') %>% mutate(sig='up')
Neutrophil.gene.up.ego.result <- Neutrophil.gene.up.ego@result %>% mutate(celltype='Neutrophil') %>% mutate(sig='up')
NK.gene.up.ego.result <- NK.gene.up.ego@result %>% mutate(celltype='NK') %>% mutate(sig='up')
Pla.gene.up.ego.result <- Pla.gene.up.ego@result %>% mutate(celltype='Pla') %>% mutate(sig='up')
T.gene.up.ego.result <- T.gene.up.ego@result %>% mutate(celltype='T') %>% mutate(sig='up')
Unciliated.gene.up.ego.result <- Unciliated.gene.up.ego@result %>% mutate(celltype='Unciliated') %>% mutate(sig='up')
SMC.gene.up.ego.result <- SMC.gene.up.ego@result %>% mutate(celltype='SMC') %>% mutate(sig='up')
Fibroblast.gene.up.ego.result <- Fibroblast.gene.up.ego@result %>% mutate(celltype='Fibroblast') %>% mutate(sig='up')



immune.upgene.ego.result <- rbind(B.gene.up.ego.result,Dendritic.gene.up.ego.result,
                                    Macrophage.gene.up.ego.result,Mast.gene.up.ego.result,
                                    Neutrophil.gene.up.ego.result,NK.gene.up.ego.result,
                                    T.gene.up.ego.result
                                    ) %>%
  filter(p.adjust < 0.05) %>% 
  dcast(Description~celltype, value.var = "p.adjust") %>% 
  mutate(number=7-rowSums(is.na(.)))

immune.downgene.ego.result <- rbind(B.gene.down.ego.result,Dendritic.gene.down.ego.result,
                                    Macrophage.gene.down.ego.result,Mast.gene.down.ego.result,
                                    Neutrophil.gene.down.ego.result,NK.gene.down.ego.result,
                                    T.gene.down.ego.result
) %>%
  filter(p.adjust < 0.05) %>% 
  dcast(Description~celltype, value.var = "p.adjust") %>% 
  mutate(number=7-rowSums(is.na(.)))



immunen.upgene.ego.result <- rbind(Ciliated.gene.up.ego.result,Endothelium.gene.up.ego.result,
                                   Pla.gene.up.ego.result,Unciliated.gene.up.ego.result,
                                   SMC.gene.up.ego.result,Fibroblast.gene.up.ego.result
) %>%
  filter(p.adjust < 0.05) %>% 
  dcast(Description~celltype, value.var = "p.adjust") %>% 
  mutate(number=7-rowSums(is.na(.)))

immunen.downgene.ego.result <- rbind(Ciliated.gene.down.ego.result,Endothelium.gene.down.ego.result,
                                   Pla.gene.down.ego.result,Unciliated.gene.down.ego.result,
                                   SMC.gene.down.ego.result,Fibroblast.gene.down.ego.result
) %>%
  filter(p.adjust < 0.05) %>% 
  dcast(Description~celltype, value.var = "p.adjust") %>% 
  mutate(number=7-rowSums(is.na(.)))

saveRDS(immune.upgene.ego.result, file='immune.upgene.ego.result.rds')
saveRDS(immune.downgene.ego.result, file='immune.downgene.ego.result.rds')
saveRDS(immunen.downgene.ego.result, file='immunen.downgene.ego.result.rds')
saveRDS(immunen.upgene.ego.result, file='immunen.upgene.ego.result.rds')

#######lgt chemokins

MTET.diffgene <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/testis/MTET/new/MTET.diffgene.rds')
## 柱状图绘制



ggplot(cellprop.melt,aes(x=celltype,y=value,fill=variable))+scale_fill_manual(values = c("lightblue","#D53E4F"))+
  stat_summary(mapping=aes(fill = variable),fun=mean,fun.args = list(mult=1),width=0.7)+
  geom_bar(stat="identity", position=position_dodge(.9))  +
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2, position=position_dodge(.9))+
  geom_jitter(aes(fill = variable),, position=position_dodge(.9),shape=21, size = 2,alpha=0.9)+
  theme_classic()



ggplot(cellprop.melt,aes(x=celltype,y=value,fill=variable))+scale_fill_manual(values = c("lightblue","#D53E4F"))+
  geom_bar(stat="identity", position=position_dodge(.9),size=1)  +
  geom_jitter(aes(fill = variable),position=position_dodge(.9),shape=21, size = 2,alpha=0.9)+
  theme_classic()

ggplot(cellprop_macrophage,aes(type,proportion))+  
  geom_bar()+ 
  theme_grey(base_size = 24)+  
  labs(x=NULL,y=NULL)+
  theme(axis.text= element_text(colour = "black",face='bold'), 
        axis.title = element_text(face = 'bold')) +
  geom_point(data=cellprop_macrophage,aes(type,proportion),size=3,pch=19)


VlnPlot(MTET.total.sct.new.Macrophage,
        features=c("Ccl8","Ccl9"),
        split.by = 'group',
        pt.size = 0.4,cols = c('tomato',"#3288BD"))

MTET.total.sct.new <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/testis/MTET/new/MTET.total.sct.new.rds')
MTET.total.sct.new.Macrophage <- subset(MTET.total.sct.new, idents='Macrophage')


apoptosis_PCR <- MTET.diffgene %>% 
  filter(gene %in% apoptosis_genesymbol$apoptosis_gene) %>% 
  filter(sig != 'non') %>% 
  filter(celltype == 'Round STids' |
           celltype == 'Elongating STids' |
           celltype == 'SPG' |
           celltype == 'Scytes' |
           celltype == 'Telocytes' |
           celltype == 'Leydig' |
           celltype == 'Macrophage' |
           celltype == 'Innate Lymph' |
           celltype == 'Sertoli'
           )

write.table(apoptosis_PCR, file='../lgt/apoptosis_PCR.txt', sep='\t',quote=F,row.names = F)



############# F
adhesion_Fibroblast.sct <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion_Fibroblast.sct.rds')

DefaultAssay(adhesion_Fibroblast.sct) <- 'RNA'
Idents(adhesion_Fibroblast.sct) <- 'integrated_snn_res.0.8'
DefaultAssay(adhesion_Fibroblast.sct) <- 'integrated'
adhesion_Fibroblast.sct <- FindClusters(adhesion_Fibroblast.sct,resolution = 0.8)

cluster_markers_all_Fibroblast <- adhesion_Fibroblast.sct.0209@misc$cluster_markers_all_Fibroblast <- Seurat::FindAllMarkers(object = adhesion_Fibroblast.sct.0209, 
                                                                                                                    assay = "RNA",
                                                                                                                    slot = "data",
                                                                                                                    verbose = TRUE, 
                                                                                                                    only.pos = TRUE, 
                                                                                                                    logfc.threshold = 0.5,
                                                                                                                    min.pct = 0.5)

cluster_markers_all_Fibroblast_40 <- cluster_markers_all_Fibroblast %>% group_by(cluster) %>% top_n(n =10, wt = avg_logFC)

adhesion_Fibroblast.sct.0209.average <- AverageExpression(adhesion_Fibroblast.sct.0209, return.seurat=TRUE )


DoHeatmap(adhesion_Fibroblast.sct.0209.average, features = unique(cluster_markers_all_Fibroblast_40$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#FFCCFF','#FF9966','#CCCCCC','#CCCCFF',
                                                  '#66CCFF','tomato','lightgreen')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 

 ### Fibroblast 注释
adhesion_Fibroblast.sct@meta.data <- adhesion_Fibroblast.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.8 == "9" ~ 'F1',
                            integrated_snn_res.0.8 == "5" ~ 'F1',
                            integrated_snn_res.0.8 == "3" ~ 'F2',
                            integrated_snn_res.0.8 == "2" ~ 'F3',
                            integrated_snn_res.0.8 == "11" ~ 'F4',
                            
                            TRUE ~ as.character(integrated_snn_res.0.8)
  ))

Idents(adhesion_Fibroblast.sct) <- 'celltype'

DefaultAssay(adhesion_Fibroblast.sct) <- 'RNA'
FeaturePlot(adhesion_Fibroblast.sct, features = c('COL1A1','KRT18','PTPRC'))


########### 剔除cluster0
adhesion_Fibroblast.sct.0209 <- subset(adhesion_Fibroblast.sct, idents=c('1', 'F2', 'F4', 'F1', 
                                                                         'F3', '6','4', '7', '8', '10'))


DefaultAssay(adhesion_Fibroblast.sct.0209) <- 'integrated'
adhesion_Fibroblast.sct.0209 <- RunPCA(adhesion_Fibroblast.sct.0209, verbose = FALSE)
ElbowPlot(adhesion_Fibroblast.sct.0209, ndims=50)
adhesion_Fibroblast.sct.0209 <- RunUMAP(adhesion_Fibroblast.sct.0209, dims = 1:50, verbose = FALSE)
adhesion_Fibroblast.sct.0209 <- RunTSNE(adhesion_Fibroblast.sct.0209, dims = 1:50, verbose = FALSE)
adhesion_Fibroblast.sct.0209 <- FindNeighbors(adhesion_Fibroblast.sct.0209, dims = 1:50)
adhesion_Fibroblast.sct.0209 <- FindClusters(adhesion_Fibroblast.sct.0209,resolution = 2)

adhesion_Fibroblast.sct.0209@meta.data <- adhesion_Fibroblast.sct.0209@meta.data %>%
  mutate(celltype=case_when(celltype == "1" ~ 'F5',
                            celltype == "8" ~ 'F5',
                            celltype == "10" ~ 'F5',
                            celltype == "4" ~ 'F6',
                            celltype == "6" ~ 'F7',
                            celltype == "7" ~ 'F8',
                            TRUE ~ as.character(celltype)
  ))

adhesion_Fibroblast.sct.0209@meta.data <- adhesion_Fibroblast.sct.0209@meta.data %>%
  mutate(celltype=case_when(celltype == "F8" ~ 'F2',
                            TRUE ~ as.character(celltype)
  ))

Idents(adhesion_Fibroblast.sct.0209) <- 'celltype'


######### fibroblast go terms

F1.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F1')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F2.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F2')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


F3.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F3')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)
 
F4.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F4')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F5.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F5')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F6.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F6')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F7.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F7')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F8.gene.ego <- cluster_markers_all_Fibroblast %>%
  dplyr::filter(cluster == 'F8')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

saveRDS(adhesion_Fibroblast.sct.0209, file='adhesion_Fibroblast.sct.0209.rds')

#########细胞比例

adhesion_Fibroblast.sct.0209@meta.data <- adhesion_Fibroblast.sct.0209@meta.data %>% 
  mutate(CellID=rownames(adhesion_Fibroblast.sct.0209@meta.data))
adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>% 
  mutate(CellID=rownames(adhesion.total.sct@meta.data))

adhesion_Fibroblast.sct.0209_annote <- adhesion_Fibroblast.sct.0209@meta.data %>% 
  dplyr::select(celltype, CellID)

adhesion.total.sct@meta.data <- left_join(adhesion.total.sct@meta.data, 
                                          adhesion_Fibroblast.sct.0209_annote, by = "CellID")

rownames(adhesion.total.sct@meta.data) <- adhesion.total.sct@meta.data$CellID


adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>%
  mutate(celltype=case_when(CellType == "FIB" ~ celltype,
                            TRUE ~ as.character(CellType)
  ))

Idents(adhesion.total.sct) <- 'celltype'

##细胞比例

###########细胞比例 整体
library(cowplot)
library(ggplot2)
library(dplyr)
library(ggpubr)

Idents(adhesion.total.sct) <- 'group'
adhesion <- subset(adhesion.total.sct, idents='adhesion')
normal <- subset(adhesion.total.sct, idents='normal')

Idents(adhesion) <- 'celltype'
Idents(normal) <- 'celltype'

cellprop <- rbind(as.data.frame(prop.table(table(Idents(adhesion), adhesion$orig.ident))),
                  as.data.frame(prop.table(table(Idents(normal), normal$orig.ident))))

colnames(cellprop)<-c("celltype","group","proportion")

cellprop$group <- as.character(cellprop$group)

cellprop <- cellprop %>% 
  mutate(group1=case_when(group == "adhesion1" ~ 'adhesion',
                          group == "adhesion2"  ~ 'adhesion',
                          group == "adhesion3" ~ 'adhesion',
                          group == "adhesion4" ~ 'adhesion',
                          group == "control1" ~ 'normal',
                          group == "control2" ~ 'normal',
                          group == "control3" ~ 'normal',
                          group == "control4" ~ 'normal',
                          TRUE ~ '1'))

ggplot(cellprop %>% filter_all(any_vars(str_detect(., pattern="Fibroblast"))), aes(celltype, weight = proportion, fill = group1)) +
  geom_bar(color = "black", width = .7, position = 'dodge') + theme_classic2() + coord_flip() + scale_fill_manual(values = c('pink','lightblue'))


saveRDS(adhesion.total.sct, file='/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion.total.sct.0210.rds')

######macrophage 注释
adhesion_Macrophage.sct <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion_Macrophage.sct.rds')

DimPlot(adhesion_Macrophage.sct,reduction = 'umap',label=TRUE,pt.size = 0.3,split.by = 'group')

adhesion_Macrophage.sct <- RunPCA(adhesion_Macrophage.sct, verbose = FALSE)
ElbowPlot(adhesion_Macrophage.sct, ndims=50)
adhesion_Macrophage.sct <- RunUMAP(adhesion_Macrophage.sct, dims = 1:40, verbose = FALSE)
adhesion_Macrophage.sct <- RunTSNE(adhesion_Macrophage.sct, dims = 1:40, verbose = FALSE)
adhesion_Macrophage.sct <- FindNeighbors(adhesion_Macrophage.sct, dims = 1:40)
adhesion_Macrophage.sct <- FindClusters(adhesion_Macrophage.sct,resolution = 0.4)

DimPlot(adhesion_Macrophage.sct.0212,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 

Idents(adhesion_Macrophage.sct) <- 'CellType'
adhesion_Macrophage.sct.0212 <- subset(adhesion_Macrophage.sct, idents='MAC')
Idents(adhesion_Macrophage.sct.0212) <- 'CellType'


adhesion_Macrophage.sct.0212 <- RunPCA(adhesion_Macrophage.sct.0212, verbose = FALSE)
ElbowPlot(adhesion_Macrophage.sct.0212, ndims=50)
adhesion_Macrophage.sct.0212 <- RunUMAP(adhesion_Macrophage.sct.0212, dims = 1:40, verbose = FALSE)
adhesion_Macrophage.sct.0212 <- RunTSNE(adhesion_Macrophage.sct.0212, dims = 1:40, verbose = FALSE)
adhesion_Macrophage.sct.0212 <- FindNeighbors(adhesion_Macrophage.sct.0212, dims = 1:40)
adhesion_Macrophage.sct.0212 <- FindClusters(adhesion_Macrophage.sct.0212,resolution = 0.4)

cluster_markers_all_Macrophage <- adhesion_Macrophage.sct.0212@misc$cluster_markers_all_Macrophage <- 
  Seurat::FindAllMarkers(object = adhesion_Macrophage.sct.0212, 
                         assay = "RNA",
                         slot = "data",
                         verbose = TRUE, 
                         only.pos = TRUE, 
                         logfc.threshold = 0.5,
                         min.pct = 0.5)

cluster_markers_all_Macrophage_40 <- cluster_markers_all_Macrophage %>% 
  group_by(cluster) %>% top_n(n =10, wt = avg_logFC)

adhesion_Macrophage.sct.0212.average <- AverageExpression(adhesion_Macrophage.sct.0212, return.seurat=TRUE )

DefaultAssay(adhesion_Macrophage.sct.0212) <- 'RNA'
DoHeatmap(adhesion_Macrophage.sct.0212.average, features = unique(cluster_markers_all_Macrophage$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#FFCCFF','#FF9966','#CCCCCC','#CCCCFF',
                                                  '#66CCFF','tomato','lightgreen')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 


####macrophage go terms


adhesion_Macrophage.sct.0212@meta.data <- adhesion_Macrophage.sct.0212@meta.data %>%
  mutate(celltype=case_when(celltype == "F0" ~ 'M0',
                            celltype == "F1" ~ 'M1',
                            celltype == "F2" ~ 'M2',
                            celltype == "F3" ~ 'M3',
                            celltype == "F4" ~ 'M4',
                            celltype == "F5" ~ 'M5',
                            celltype == "F6" ~ 'M6',
                            celltype == "F7" ~ 'M7',
                            TRUE ~ as.character(CellType)
  ))






Idents(adhesion_Macrophage.sct.0212) <- 'celltype'

######### Macrophage go terms

F1.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F1')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F2.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F2')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


F3.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F3')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F4.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F4')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F5.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F5')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F6.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F6')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F7.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F7')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

F0.gene.ego <- cluster_markers_all_Macrophage %>%
  dplyr::filter(cluster == 'F0')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


##############
adhesion_Macrophage.sct.0212@meta.data <- adhesion_Macrophage.sct.0212@meta.data %>% 
  mutate(CellID=rownames(adhesion_Macrophage.sct.0212@meta.data))

adhesion.total.sct@meta.data <- adhesion.total.sct@meta.data %>% 
  mutate(CellID=rownames(adhesion.total.sct@meta.data))

adhesion_Macrophage.sct.0212_annote <- adhesion_Macrophage.sct.0212@meta.data %>% 
  dplyr::select(celltype, CellID)

adhesion.total.sct.0210@meta.data <- left_join(adhesion.total.sct.0210@meta.data, 
                                          adhesion_Macrophage.sct.0212_annote, by = "CellID")

rownames(adhesion.total.sct.0210@meta.data) <- adhesion.total.sct.0210@meta.data$CellID


adhesion.total.sct.0210@meta.data <- adhesion.total.sct.0210@meta.data %>%
  mutate(celltype=case_when(CellType == "MAC" ~ celltype.y,
                            TRUE ~ as.character(celltype.x)
  ))

Idents(adhesion.total.sct.0210) <- 'celltype'

DimPlot(adhesion.total.sct.0210,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 

adhesion.total.sct.0210 <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion.total.sct.0210.rds')

rownames(adhesion.total.sct.0210@meta.data) <- adhesion.total.sct.0210@meta.data$CellID


########### Nk cell

adhesion_NK <- subset(adhesion.total.sct.0210,idents = c('NK'))

########Macrophage细胞注释：

adhesion_NK.list <- SplitObject(adhesion_NK, split.by = "orig.ident")

adhesion_NK.list <- lapply(X = adhesion_NK.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = adhesion_NK.list)
adhesion_NK.list <- PrepSCTIntegration(object.list = adhesion_NK.list, anchor.features = features)
adhesion_NK.list.anchors <- FindIntegrationAnchors(object.list = adhesion_NK.list, normalization.method = "SCT", 
                                                  anchor.features = features)
adhesion_NK.sct <- IntegrateData(anchorset = adhesion_NK.list.anchors, normalization.method = "SCT")


adhesion_NK.sct <- RunPCA(adhesion_NK.sct, verbose = FALSE)
ElbowPlot(adhesion_NK.sct, ndims=50)
adhesion_NK.sct <- RunUMAP(adhesion_NK.sct, dims = 1:40, verbose = FALSE)
adhesion_NK.sct <- RunTSNE(adhesion_NK.sct, dims = 1:40, verbose = FALSE)
adhesion_NK.sct <- FindNeighbors(adhesion_NK.sct, dims = 1:40)
adhesion_NK.sct <- FindClusters(adhesion_NK.sct,resolution = 0.8)


DimPlot(adhesion_NK.sct,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 


cluster_markers_all_NK <- adhesion_NK.sct@misc$cluster_markers_all_NK <- 
  Seurat::FindAllMarkers(object = adhesion_NK.sct, 
                         assay = "RNA",
                         slot = "data",
                         verbose = TRUE, 
                         only.pos = TRUE, 
                         logfc.threshold = 0.5,
                         min.pct = 0.5)

cluster_markers_all_NK_40 <- cluster_markers_all_NK %>% 
  group_by(cluster) %>% top_n(n =10, wt = avg_logFC)

adhesion_NK.sct.average <- AverageExpression(adhesion_NK.sct, return.seurat=TRUE )

DefaultAssay(adhesion_NK.sct) <- 'RNA'
DefaultAssay(adhesion_NK.sct) <- 'integrated'

DoHeatmap(adhesion_NK.sct.average, features = unique(cluster_markers_all_NK$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#FFCCFF','#FF9966','#CCCCCC','#CCCCFF',
                                                  '#66CCFF','tomato','lightgreen')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 


adhesion_NK.sct@meta.data <- adhesion_NK.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.4 == "0" ~ 'N0',
                            integrated_snn_res.0.4 == "1" ~ 'N1',
                            integrated_snn_res.0.4 == "2" ~ 'N2',
                            integrated_snn_res.0.4 == "3" ~ 'N2',
                            integrated_snn_res.0.4 == "4" ~ 'N1',
                            integrated_snn_res.0.4 == "5" ~ 'N3',
                            TRUE ~ as.character(CellType)
  ))

adhesion_NK.sct@meta.data <- adhesion_NK.sct@meta.data %>%
  mutate(celltype=case_when(celltype == "N0" ~ 'N0',
                            celltype == "N2" ~ 'N0',
                            TRUE ~ as.character(celltype)
  ))


Idents(adhesion_NK.sct) <- 'celltype'


adhesion_NK.sct@meta.data <- adhesion_NK.sct@meta.data %>% 
  mutate(CellID=rownames(adhesion_NK.sct@meta.data))

adhesion.total.sct.0210@meta.data <- adhesion.total.sct.0210@meta.data %>% 
  mutate(CellID=rownames(adhesion.total.sct@meta.data))

adhesion_NK.sct_annote <- adhesion_NK.sct@meta.data %>% 
  dplyr::select(celltype, CellID)

adhesion.total.sct.0210@meta.data <- left_join(adhesion.total.sct.0210@meta.data, 
                                               adhesion_NK.sct_annote, by = "CellID")

rownames(adhesion.total.sct.0210@meta.data) <- adhesion.total.sct.0210@meta.data$CellID


adhesion.total.sct.0210@meta.data <- adhesion.total.sct.0210@meta.data %>%
  mutate(celltype=case_when(CellType == "NK" ~ celltype.y.y,
                            TRUE ~ as.character(celltype.x.x)
  ))

Idents(adhesion.total.sct.0210) <- 'celltype'

########### T cell analys
adhesion_T.sct <- readRDS('/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion_T.sct.rds')

DimPlot(adhesion_T.sct,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') 

DefaultAssay(adhesion_T.sct) <- 'RNA'

DotPlot(adhesion_T.sct, features = c('CD4','CD8A'), 
        cols =  c("white", "deeppink3"),scale.by='size') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_text(size=12,family="Arial")) +
  theme(axis.text.x=element_text(size=12,angle=90,family="Arial"),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.title.x=element_text(size = 18,family="Arial"),
        axis.title.y=element_text(size = 18,family="Arial")) + ylab('Celltype') + xlab('Marker gene') +
  theme(legend.key.height = unit(10,'pt'))

FeaturePlot(adhesion_T.sct, features = c('CD4','CD8A',
                                         'CD3D','NKG7','TRDC'))

DotPlot(adhesion_T.sct, features = c('CTLA4','CCR7','LEF1','TCF7','CXCR6','CD40LG','TRDC'), 
        cols =  c("white", "deeppink3"),scale.by='size') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_text(size=12,family="Arial")) +
  theme(axis.text.x=element_text(size=12,angle=90,family="Arial"),
        axis.text.y=element_text(size=12,family="Arial"),
        axis.title.x=element_text(size = 18,family="Arial"),
        axis.title.y=element_text(size = 18,family="Arial")) + ylab('Celltype') + xlab('Marker gene') +
  theme(legend.key.height = unit(10,'pt'))

DefaultAssay(adhesion_T.sct) <- 'integrated'
adhesion_T.sct <- RunPCA(adhesion_T.sct, verbose = FALSE)
ElbowPlot(adhesion_T.sct, ndims=50)
adhesion_T.sct <- RunUMAP(adhesion_T.sct, dims = 1:50, verbose = FALSE)
adhesion_T.sct <- RunTSNE(adhesion_T.sct, dims = 1:50, verbose = FALSE)
adhesion_T.sct <- FindNeighbors(adhesion_T.sct, dims = 1:50)
adhesion_T.sct <- FindClusters(adhesion_T.sct,resolution = 0.4)

cluster_markers_all_T <- adhesion_T.sct@misc$cluster_markers_all_T <- 
  Seurat::FindAllMarkers(object = adhesion_T.sct, 
                         assay = "RNA",
                         slot = "data",
                         verbose = TRUE, 
                         only.pos = TRUE, 
                         logfc.threshold = 0.5,
                         min.pct = 0.5)

cluster_markers_all_T_40 <- cluster_markers_all_T %>% 
  group_by(cluster) %>% top_n(n =10, wt = avg_logFC)

adhesion_T.sct.average <- AverageExpression(adhesion_T.sct, return.seurat=TRUE )

DefaultAssay(adhesion_T.sct) <- 'RNA'
DoHeatmap(adhesion_T.sct.average, features = unique(cluster_markers_all_T$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#FFCCFF','#FF9966','#CCCCCC','#CCCCFF',
                                                  '#66CCFF','tomato','lightgreen')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 


adhesion_T.sct@meta.data <- adhesion_T.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.4 == "1" ~ 'CD4',
                            integrated_snn_res.0.4 == "7" ~ 'CD4',
                            integrated_snn_res.0.4 == "8" ~ 'CD4',
                            integrated_snn_res.0.4 == "0" ~ 'CD8',
                            integrated_snn_res.0.4 == "2" ~ 'CD8',
                            integrated_snn_res.0.4 == "3" ~ 'CD8',
                            integrated_snn_res.0.4 == "4" ~ 'CD8',
                            integrated_snn_res.0.4 == "5" ~ 'CD8',
                            integrated_snn_res.0.4 == "11" ~ 'CD8',
                            integrated_snn_res.0.4 == "12" ~ 'CD8',
                            integrated_snn_res.0.4 == "10" ~ 'CD8',
                            integrated_snn_res.0.4 == "6" ~ 'NKT',
                            integrated_snn_res.0.4 == "9" ~ 'γδT',
                            TRUE ~ as.character(integrated_snn_res.0.4)
  ))
Idents(adhesion_T.sct) <- 'celltype'

DimPlot(adhesion_T.sct,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 


adhesion_T.sct@meta.data <- adhesion_T.sct@meta.data %>% 
  mutate(CellID=rownames(adhesion_T.sct@meta.data))

adhesion_T.sct_annote <- adhesion_T.sct@meta.data %>% 
  dplyr::select(celltype, CellID)

adhesion.total.sct.0210@meta.data <- left_join(adhesion.total.sct.0210@meta.data, 
                                               adhesion_T.sct_annote, by = "CellID")

rownames(adhesion.total.sct.0210@meta.data) <- adhesion.total.sct.0210@meta.data$CellID


adhesion.total.sct.0210@meta.data <- adhesion.total.sct.0210@meta.data %>%
  mutate(celltype=case_when(CellType == 'T' ~ celltype.y.y.y,
                            TRUE ~ as.character(celltype.x.x.x)
  ))

Idents(adhesion.total.sct.0210) <- 'celltype'


########### Nk cell

Epi.sct <- subset(adhesion.total.sct.0210,idents = c('UNCILIA','CILIA'))

DefaultAssay(Epi.sct) <- 'integrated'
Epi.sct <- RunPCA(Epi.sct, verbose = FALSE)
ElbowPlot(Epi.sct, ndims=50)
Epi.sct <- RunUMAP(Epi.sct, dims = 1:40, verbose = FALSE)
Epi.sct <- RunTSNE(Epi.sct, dims = 1:40, verbose = FALSE)
Epi.sct <- FindNeighbors(Epi.sct, dims = 1:40)
Epi.sct <- FindClusters(Epi.sct,resolution = 0.4)

DimPlot(Epi.sct,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 


cluster_markers_all_Epi <- Epi.sct@misc$cluster_markers_all_Epi <- 
  Seurat::FindAllMarkers(object = Epi.sct, 
                         assay = "RNA",
                         slot = "data",
                         verbose = TRUE, 
                         only.pos = TRUE, 
                         logfc.threshold = 0.5,
                         min.pct = 0.5)

cluster_markers_all_T_40 <- cluster_markers_all_T %>% 
  group_by(cluster) %>% top_n(n =10, wt = avg_logFC)

Epi.sct.average <- AverageExpression(Epi.sct, return.seurat=TRUE )

DefaultAssay(Epi.sct) <- 'RNA'
DoHeatmap(Epi.sct.average, features = unique(cluster_markers_all_Epi$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#FFCCFF','#FF9966','#CCCCCC','#CCCCFF',
                                                  '#66CCFF','tomato','lightgreen')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 


Idents(Epi.sct) <- 'integrated_snn_res.0.4'


UNCILIA.sct <- subset(adhesion.total.sct.0210,idents = c('UNCILIA'))
CILIA.sct <- subset(adhesion.total.sct.0210,idents = c('CILIA'))

DefaultAssay(CILIA.sct) <- 'integrated'
CILIA.sct <- RunPCA(CILIA.sct, verbose = FALSE)
ElbowPlot(Epi.sct, ndims=50)
CILIA.sct <- RunUMAP(CILIA.sct, dims = 1:50, verbose = FALSE)
CILIA.sct <- RunTSNE(CILIA.sct, dims = 1:50, verbose = FALSE)
CILIA.sct <- FindNeighbors(CILIA.sct, dims = 1:50)
CILIA.sct <- FindClusters(CILIA.sct,resolution = 0.4)

DefaultAssay(UNCILIA.sct) <- 'integrated'
UNCILIA.sct <- RunPCA(UNCILIA.sct, verbose = FALSE)
ElbowPlot(Epi.sct, ndims=50)
UNCILIA.sct <- RunUMAP(UNCILIA.sct, dims = 1:50, verbose = FALSE)
UNCILIA.sct <- RunTSNE(UNCILIA.sct, dims = 1:50, verbose = FALSE)
UNCILIA.sct <- FindNeighbors(UNCILIA.sct, dims = 1:50)
UNCILIA.sct <- FindClusters(UNCILIA.sct,resolution = 0.4)

DimPlot(UNCILIA.sct,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 

cluster_markers_all_UNCILIA <- UNCILIA.sct@misc$cluster_markers_all_UNCILIA <- 
  Seurat::FindAllMarkers(object = UNCILIA.sct, 
                         assay = "RNA",
                         slot = "data",
                         verbose = TRUE, 
                         only.pos = TRUE, 
                         logfc.threshold = 0.5,
                         min.pct = 0.5)

cluster_markers_all_UNCILIA_40 <- cluster_markers_all_UNCILIA %>% 
  group_by(cluster) %>% top_n(n =10, wt = avg_logFC)

UNCILIA.sct.average <- AverageExpression(UNCILIA.sct, return.seurat=TRUE )

DefaultAssay(UNCILIA.sct) <- 'RNA'
DoHeatmap(UNCILIA.sct.average, features = unique(cluster_markers_all_UNCILIA$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#FFCCFF','#FF9966','#CCCCCC','#CCCCFF',
                                                  '#66CCFF','tomato','lightgreen')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 


UNCILIA.sct@meta.data <- UNCILIA.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.4 == "0" ~ 'U1',
                            integrated_snn_res.0.4 == "10" ~ 'U1',
                            integrated_snn_res.0.4 == "9" ~ 'U1',
                            TRUE ~ as.character(integrated_snn_res.0.4)
  ))

Idents(UNCILIA.sct) <- 'celltype'

####### go terms

U1.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == 'U1')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

U2.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '2')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

U3.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '3')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

U7.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '7')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

U5.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '5')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


U6.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '6')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

U8.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '8')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

U4.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '4')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

U0.gene.ego <- cluster_markers_all_UNCILIA %>%
  dplyr::filter(cluster == '1')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


DimPlot(CILIA.sct,reduction = 'tsne',label=TRUE,pt.size = 0.3,split.by = 'group') + 
  scale_colour_manual(values =  c('lightblue','#FFCCFF','#FF9966','pink','#CCFF00','#CCCCFF',
                                  '#CCCCCC','mediumpurple','#66FF00','#EE0000','#66CCFF','tomato'))+
  guides(color=guide_legend(override.aes = list(size=6))) 



DefaultAssay(CILIA.sct) <- 'RNA'
cluster_markers_all_CILIA <- CILIA.sct@misc$cluster_markers_all_CILIA <- 
  Seurat::FindAllMarkers(object = CILIA.sct, 
                         assay = "RNA",
                         slot = "data",
                         verbose = TRUE, 
                         only.pos = TRUE, 
                         logfc.threshold = 0.5,
                         min.pct = 0.5)

cluster_markers_all_UNCILIA_40 <- cluster_markers_all_UNCILIA %>% 
  group_by(cluster) %>% top_n(n =10, wt = avg_logFC)

CILIA.sct.average <- AverageExpression(CILIA.sct, return.seurat=TRUE )


DoHeatmap(CILIA.sct.average, features = unique(cluster_markers_all_CILIA$gene), angle = 45,
          size = 4, draw.lines = F,group.colors=c('#FFCCFF','#FF9966','#CCCCCC','#CCCCFF',
                                                  '#66CCFF','tomato','lightgreen')) + 
  scale_fill_gradientn(colors = c('lightblue', 'white', 'deeppink')) +  
  theme(legend.text = element_text(colour = 'black', size = 12,)) +
  guides(color=guide_legend(override.aes = list(size=8))) 


CILIA.sct@meta.data <- CILIA.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.4 == "0" ~ 'C1',
                            integrated_snn_res.0.4 == "1" ~ 'C1',
                            integrated_snn_res.0.4 == "2" ~ 'C2',
                            integrated_snn_res.0.4 == "3" ~ 'C3',
                            TRUE ~ as.character(CellType)
  ))

Idents(CILIA.sct) <- 'celltype'


#### cell annotes

C1.gene.ego <- cluster_markers_all_CILIA %>%
  dplyr::filter(cluster == 'C1')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

C2.gene.ego <- cluster_markers_all_CILIA %>%
  dplyr::filter(cluster == 'C2')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

C3.gene.ego <- cluster_markers_all_CILIA %>%
  dplyr::filter(cluster == 'C3')  %>%
  dplyr::pull(gene, name=avg_logFC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


UNCILIA.sct@meta.data <- UNCILIA.sct@meta.data %>%
  mutate(celltype=case_when(integrated_snn_res.0.4 == "0" ~ 'U1',
                            integrated_snn_res.0.4 == "10" ~ 'U1',
                            integrated_snn_res.0.4 == "9" ~ 'U1',
                            integrated_snn_res.0.4 == "1" ~ 'U0',
                            integrated_snn_res.0.4 == "2" ~ 'U2',
                            integrated_snn_res.0.4 == "3" ~ 'U3',
                            integrated_snn_res.0.4 == "4" ~ 'U4',
                            integrated_snn_res.0.4 == "5" ~ 'U5',
                            integrated_snn_res.0.4 == "6" ~ 'U6',
                            integrated_snn_res.0.4 == "7" ~ 'U7',
                            integrated_snn_res.0.4 == "8" ~ 'U8',
                            TRUE ~ as.character(integrated_snn_res.0.4)
  ))

Idents(UNCILIA.sct) <- 'celltype'

UNCILIA.sct@meta.data <- UNCILIA.sct@meta.data %>% 
  mutate(CellID=rownames(UNCILIA.sct@meta.data))

UNCILIA.sct_annote <- UNCILIA.sct@meta.data %>% 
  dplyr::select(celltype, CellID)

adhesion.total.sct.0210@meta.data <- left_join(adhesion.total.sct.0210@meta.data, 
                                               UNCILIA.sct_annote, by = "CellID")

rownames(adhesion.total.sct.0210@meta.data) <- adhesion.total.sct.0210@meta.data$CellID


adhesion.total.sct.0210@meta.data <- adhesion.total.sct.0210@meta.data %>%
  mutate(celltype=case_when(CellType == 'UNCILIA' ~ celltype.y.y.y.y,
                            TRUE ~ as.character(celltype.x.x.x.x)
  ))

CILIA.sct@meta.data <- CILIA.sct@meta.data %>% 
  mutate(CellID=rownames(CILIA.sct@meta.data))

CILIA.sct_annote <- CILIA.sct@meta.data %>% 
  dplyr::select(celltype, CellID)

adhesion.total.sct.0210@meta.data <- left_join(adhesion.total.sct.0210@meta.data, 
                                               CILIA.sct_annote, by = "CellID")

rownames(adhesion.total.sct.0210@meta.data) <- adhesion.total.sct.0210@meta.data$CellID


adhesion.total.sct.0210@meta.data <- adhesion.total.sct.0210@meta.data %>%
  mutate(celltype=case_when(CellType == 'CILIA' ~ celltype.y.y.y.y.y,
                            TRUE ~ as.character(celltype.x.x.x.x.x)
  ))
Idents(adhesion.total.sct.0210) <- 'celltype'



adhesion.total.sct.new <- subset(adhesion.total.sct.0210, 
                                 idents = c("CD8","U2","U0","M0","N0",
                                            "C1", "F5","MONO","U7",  
                                            "F2", "U1", "F4", "U6",  
                                            "U8", "U4","MAST","C3", 
                                            "B","ENDO", "U3", "N1",
                                            "M2","M1","N3","U5","F1",   
                                            "CD4","F3","γδT","F7", "M4","F6",   
                                            "F8", "C2","NKT","M3", "DC","M5",   
                                           "M6","M7"))

Idents(adhesion.total.sct.new) <- 'celltype'


DefaultAssay(adhesion.total.sct.new) <- 'integrated'
adhesion.total.sct.new <- RunPCA(adhesion.total.sct.new, verbose = FALSE)
ElbowPlot(adhesion.total.sct.new, ndims=50)
adhesion.total.sct.new <- RunUMAP(adhesion.total.sct.new, dims = 1:50, verbose = FALSE)
adhesion.total.sct.new <- RunTSNE(adhesion.total.sct.new, dims = 1:50, verbose = FALSE)
adhesion.total.sct.new <- FindNeighbors(adhesion.total.sct.new, dims = 1:50)
adhesion.total.sct.new <- FindClusters(adhesion.total.sct.new,resolution = 0.4)

saveRDS(adhesion.total.sct.new, 
        file='/szrmyy/wangjgLab/scRNA/xiasy/adhesion/adhesion.total.sct.new.rds')
