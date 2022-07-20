library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(monocle3)
library(SeuratWrappers)

##load data day0_2 merge
load("./data/scMuscle_pref_seurat_day0_2_f1_6K_filtered.RData")

##figure1 m
DimPlot(object = scMuscle.pref.seurat, reduction = "umap_harmony", pt.size = .1,group.by ="RNA_res.0.25")+labs(title = "Day0+Day2")

###pseudo time
scMuscle.pref.seurat@reductions$UMAP<-scMuscle.pref.seurat@reductions$umap_harmony
hamy_day0_2_sample_cds<-as.cell_data_set(scMuscle.pref.seurat)
hamy_day0_2_sample_cds<- cluster_cells(hamy_day0_2_sample_cds,reduction = "UMAP",k = 100)
hamy_day0_2_sample_cds <- learn_graph(hamy_day0_2_sample_cds, close_loop = F,use_partition = T,learn_graph_control =list(minimal_branch_len=5))

plot_cells(hamy_day0_2_sample_cds, label_groups_by_cluster = T, label_leaves = F, label_branch_points = T,graph_label_size = 3)

##figure1 n
library(grid)
library(gridExtra)
library(ArchR)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

col.ls <- ArchRPalettes[16:30]

FeaturePlot(scMuscle.pref.seurat,features = c("Pou3f1",
                                              "Sox4",
                                              "Tbx3",
                                              "Nanog",
                                              "Sox2"), reduction = 'umap_harmony')&scale_colour_gradientn(colours =col.ls[8]$comet)

##figure1 o
FeaturePlot(scMuscle.pref.seurat,features = c("Gata4",
                                              "Hnf1b"), reduction = 'umap_harmony')&scale_colour_gradientn(colours =col.ls[8]$comet)

##sfigure1 F G I

FeaturePlot(scMuscle.pref.seurat,features = c("Ulk1",
                                              "Gabarapl2",
                                              "Otx2",
                                              "Sall2",
                                              "Gata6",
                                              "Foxa2"), reduction = 'umap_harmony')&scale_colour_gradientn(colours =col.ls[8]$comet)

## figure1 h
load(file="./data/agg_day0_day2/comb_day0_day2_final.RData")
library(grid)
library(gridExtra)

DimPlot(object = integrated, label = F)


##figure1 p
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(patchwork)
set.seed(1234)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library("BSgenome.Mmusculus.UCSC.mm10")

load("./data/comb_day0_day2_final_tbx3s.RData")
integrated_tbx3s<-integrated
load("./data/comb_day0_day2_final.RData")

###Overall "Gata4"
bg_peaks <- FindAllMarkers(
  object = integrated,
  only.pos = TRUE,
  test.use = 'LR',
  logfc.threshold = 0.15,
  min.pct = 0.1,
  latent.vars = 'nCount_peaks'
)

da_peaks <- FindMarkers(
  object = integrated,
  ident.1 = c("3"),
  ident.2 = c("0","1","2"),
  only.pos = TRUE,
  test.use = 'LR',
  logfc.threshold = 0.15,
  min.pct = 0.1,
  latent.vars = 'nCount_peaks'
)

enriched.motifs <- FindMotifs(
  object = integrated_tbx3s,
  features = rownames(da_peaks),
  background =rownames(integrated)
)

PlotFootprint(integrated, features = c("GATA4"))

MotifPlot(
  object = integrated,
  motifs = "MA0482.2"
)

###Overall "Hnf1b"
PlotFootprint(integrated, features = c("HNF1B"))
MotifPlot(
  object = integrated,
  motifs = "MA0153.2"
)


##Sfigure 1 I


tiff("Gata6_all_motif.tiff", units="in", width=10, height=5, res=300)
MotifPlot(
  object = integrated,
  motifs = "MA1104.2"
)
dev.off()

tiff("Foxa2_all_motif.tiff", units="in", width=10, height=5, res=300)
MotifPlot(
  object = integrated,
  motifs = "MA0047.3"
)
dev.off()

###Footprint on Specific motif
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
motif.matrix <- CreateMotifMatrix(
  features = granges(integrated),
  pwm = pfm,
  genome = 'mm10'
)
HNF1B_coord<-motif.matrix[motif.matrix[,"MA0153.2"]==1,"MA0153.2"]
integrated_HNF1B <- subset(x = integrated, features=names(HNF1B_coord))
integrated_HNF1B <- AddMotifs(integrated_HNF1B, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
integrated_HNF1B  <- Footprint(
  object = integrated_HNF1B ,
  motif.name = c("GATA4"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)
enriched.motifs_HNF1B <- FindMotifs(
  object = integrated_HNF1B,
  features = rownames(da_peaks_HNF1B),
  background =rownames(integrated_HNF1B)
)
PlotFootprint(integrated_HNF1B, features = c("GATA4"))

tiff("Gata4_tss_motif.tiff", units="in", width=10, height=5, res=300)
MotifPlot(
  object = integrated_HNF1B,
  motifs = "MA0482.2"
)
dev.off()

cata4_coord<-motif.matrix[motif.matrix[,"MA0482.2"]==1,"MA0482.2"]
integrated_gata4 <-integrated[names(cata4_coord),]
integrated_gata4 <- AddMotifs(integrated_gata4, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
integrated_gata4  <- Footprint(
  object = integrated_gata4 ,
  motif.name = c("HNF1B"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)

da_peaks_gata4 <- FindMarkers(
  object = integrated_gata4,
  ident.1 = c("3"),
  ident.2 = c("0","1","2"),
  only.pos = TRUE,
  test.use = 'LR',
  logfc.threshold = 0.15,
  min.pct = 0.1,
  latent.vars = 'nCount_peaks'
)
enriched.motifs_gata4 <- FindMotifs(
  object = integrated_gata4,
  features = rownames(da_peaks_gata4),
  background =rownames(integrated_gata4)
)
PlotFootprint(integrated_gata4, features = c("HNF1B"))
tiff("Hnf1b_tss_motif.tiff", units="in", width=10, height=5, res=300)
MotifPlot(
  object = integrated_gata4,
  motifs = "MA0153.2"
)
dev.off()
