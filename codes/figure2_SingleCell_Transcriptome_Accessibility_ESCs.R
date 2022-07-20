## figure2 c
library(dplyr)
library(Seurat)

# for plotting
library(ggplot2)
library(patchwork)
set.seed(1234)
library(monocle3)
library(SeuratWrappers)

load(file="./data/day6_filtered_seur_f2k_rd_cellcycle.RData")

DimPlot(object = day6, reduction = 'umap',label=F)+labs(title = "day6")

### pseudo time
DefaultAssay(day6) <- "RNA"
day6_cds<-as.cell_data_set(day6)
day6_cds<- cluster_cells(day6_cds,reduction = "UMAP",k = 30,resolution = 0.00012)
day6_cds <- learn_graph(day6_cds, close_loop = F,use_partition = T,learn_graph_control =list(minimal_branch_len=5))
plot_cells(day6_cds, label_groups_by_cluster = T, label_leaves = F, label_branch_points = T,graph_label_size = 3)
day6.min.umap <- which.min(unlist(FetchData(day6, "UMAP_2")))
day6.min.umap <- colnames(day6)[day6.min.umap]
day6_cds <- order_cells(day6_cds, root_cells = day6.min.umap)
day6_cds <- order_cells(day6_cds)
plot_cells(day6_cds, color_cells_by = "pseudotime", label_cell_groups =T, label_leaves = F, 
           label_branch_points = F,show_trajectory_graph = T,graph_label_size = 3,label_groups_by_cluster = T)

## figure2 d-h
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(ggplot2)
set.seed(1234)

library(ArchR)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

load(file="./data/mamd_day6_5k/day6_atac_5k__new_final.RData")
day6_atac_5k <- AddMetaData(object = day6_atac_5k, metadata = predicted.labels)

pPSM<-read.csv("pPSM_f.txt",header =F)[,1]
aPSM<-read.csv("aPSM_f.txt",header =F)[,1]
glyco<-read.csv("glyco.txt",header =F)[,1]
blood<-read.csv("blood.txt",header =F)[,1]
oxi<-read.csv("phos.txt",header =F)[,1]
Neuronal<-c("Nes",
            "Hes5",
            "Apoe",
            "Sox11",
            "Sox2",
            "Pax6")

del_list<-c('Rora','Tbx18','Vtn','Sox6','Ripply2','Pgm5','Psck5','Myocd','Met','Fhf18','Dmrt2','Cer1','Foxp1',
            'Zic5','Sim1','Pcsk5','Pappa','Meox2','Gadd45b','Fgf18','Cd36','Rock2','Rock1','Jup','Fyn','Calm1',
            'Actb','Apoe')
names(pPSM)<-pPSM
names(aPSM)<-aPSM
names(glyco)<-glyco
names(blood)<-blood
names(oxi)<-oxi
names(Neuronal)<-Neuronal



pPSM<-pPSM[setdiff(names(pPSM),del_list )]
aPSM<-aPSM[setdiff(names(aPSM),del_list )]
glyco<-glyco[setdiff(names(glyco),del_list )]
blood<-blood[setdiff(names(blood),del_list )]
oxi<-oxi[setdiff(names(oxi),del_list )]
Neuronal<-Neuronal[setdiff(names(Neuronal),del_list )]
glyco<-setdiff(glyco,Neuronal)

pPSM<-unname(pPSM)
aPSM<-unname(aPSM)
glyco<-unname(glyco)
blood<-unname(blood)
oxi<-unname(oxi)
Neuronal<-unname(Neuronal)

DefaultAssay(day6_atac_5k)<-"RNA"
Idents(day6_atac_5k) <- factor(day6_atac_5k$predicted.id,levels = c("0","1","2","3"))

#posterior
p1<-DotPlot(day6_atac_5k, features = pPSM, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0 ,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average estimated expression(z-score)'))

s1<-DotPlot(day6, features = pPSM, cols = c("blue", "red"), dot.scale = 8, assay="RNA",scale.min =0, scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average expression(z-score)'))


p2<-DotPlot(day6_atac_5k, features = aPSM, cols = c("blue", "red"), dot.scale = 8, assay="RNA",scale.min =0, scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average estimated expression(z-score)'))
s2<-DotPlot(day6, features = aPSM, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average expression(z-score)'))


p3<-DotPlot(day6_atac_5k, features = glyco, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average estimated expression(z-score)'))

s3<-DotPlot(day6, features = glyco, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average expression(z-score)'))

p4<-DotPlot(day6_atac_5k, features = blood, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average estimated expression(z-score)'))

s4<-DotPlot(day6, features = blood, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average expression(z-score)'))

p5<-DotPlot(day6_atac_5k, features = oxi, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average estimated expression(z-score)'))

s5<-DotPlot(day6, features = oxi, cols = c("blue", "red"), dot.scale = 8, assay="RNA",scale.min =0, scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average expression(z-score)'))


p6<-DotPlot(day6_atac_5k, features = Neuronal, cols = c("blue", "red"), dot.scale = 8, assay="RNA",scale.min =0, scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average estimated expression(z-score)'))

s6<-DotPlot(day6, features = Neuronal, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average expression(z-score)'))

library(grid)
library(gridExtra)

tiff("posterior.tiff", units="in", width=15, height=10, res=300)
grid.arrange(s1, p1, ncol = 2,widths = c(8,9))
dev.off()
tiff("Anterior.tiff", units="in", width=15, height=10, res=300)
grid.arrange(s2, p2, ncol = 2,widths = c(8,9))
dev.off()
tiff("Glycose.tiff", units="in", width=15, height=10, res=300)
grid.arrange(s3, p3, ncol = 2,widths = c(8,9),heights=c(6,6))
dev.off()
tiff("Blood.tiff", units="in", width=15, height=10, res=300)
grid.arrange(s4, p4, ncol = 2,widths = c(8,9),heights=c(6,6))
dev.off()
tiff("Oxiphose.tiff", units="in", width=15, height=10, res=300)
grid.arrange(s5, p5, ncol = 2,widths = c(8,9),heights=c(6,6))
dev.off()
tiff("Neuronal.tiff", units="in", width=15, height=10, res=300)
grid.arrange(s6, p6, ncol = 2,widths = c(8,9),heights=c(6,6))
dev.off()


## figure2 i
load("./data/scMuscle_pref_seurat_day0_2_f1_6K_filtered.RData")
day0_sep<-subset(x = scMuscle.pref.seurat, subset = sample=="Day0")

library(ggplot2)
VlnPlot(
  object = day0_sep,
  features = 'Sox2',
  pt.size = 0.1
)+ylab( "Average Expression")+geom_boxplot(width=0.1, fill="white")

###aPSM
VlnPlot(
  object = day6,
  features = 'Sox2',
  y.max=3,
  pt.size = 0.1
)+ylab( "Average Expression")+geom_boxplot(width=0.1, fill="white")

##figure2 j
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(patchwork)
set.seed(1234)

load("./data/agg_day0_day6/comb_day0_day6_revision.RData")

library(ggplot2)
DimPlot(object = integrated, label = F)

source("custom_seurat_functions.R")
plot_integrated_clusters(integrated) 

##figure2 k
tiff("sox2_day0_6_track.tiff", units="in", width=10, height=5, res=300)
CoveragePlot(
  object =integrated,
  region = "Sox2",
  features = "Sox2",
  expression.assay = "RNA",
  extend.upstream = 100,
  extend.downstream = 120000,
  links = F
)
dev.off()
tiff("sox2_day0_6_expr.tiff", units="in", width=5, height=10, res=300)
ExpressionPlot(
  object = integrated,
  features = "Sox2",
  group.by = "dataset",
  assay = "RNA"
)
dev.off()


##sfigure2
genes<-c("Pax6","Sox2","Pou3f2","Tbr1","Pgk1","Pkm","Ldha","Eno1","Fgf8")
FeaturePlot(day6,features = genes[1:3], min.cutoff = "q10",max.cutoff = "q90")&scale_colour_gradientn(colours =col.ls[8]$comet)
FeaturePlot(day6_atac_5k,features = genes[1:3], min.cutoff = "q10",max.cutoff = "q90")&scale_colour_gradientn(colours =col.ls[8]$comet)
FeaturePlot(day6,features = genes[4:6], min.cutoff = "q10",max.cutoff = "q90")&scale_colour_gradientn(colours =col.ls[8]$comet)
FeaturePlot(day6_atac_5k,features = genes[4:6], min.cutoff = "q10",max.cutoff = "q90")&scale_colour_gradientn(colours =col.ls[8]$comet)
FeaturePlot(day6,features = genes[7:9], min.cutoff = "q10",max.cutoff = "q90")&scale_colour_gradientn(colours =col.ls[8]$comet)
FeaturePlot(day6_atac_5k,features = genes[7:9], min.cutoff = "q10",max.cutoff = "q90")&scale_colour_gradientn(colours =col.ls[8]$comet)

###sfigure2 E
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(patchwork)
set.seed(1234)

load("./data/agg_day0_day6/comb_day0_day6_revision.RData")

DefaultAssay(integrated) <- 'RNA'
avg.integrated<- AverageExpression(integrated)
avg.integrated_peaks<-log1p(as.data.frame(avg.integrated$peaks))
avg.integrated_RNA<-log1p(as.data.frame(avg.integrated$RNA))
sampleinfo<-data.frame(status=c("Naïve","aPSM"))
library(ggfortify)
pcDat <- prcomp(t(avg.integrated_RNA))
# plot PCA
autoplot(pcDat,data=sampleinfo,colour="status")
