library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(patchwork)
set.seed(1234)

##figure4 b
load(file="./data/agg_day6_day8/day6_8_atac.RData")

DimPlot(object = day6_8_atac, label = F,reduction = 'umap') +labs(title = "day6+day8 scATAC")
###separate
cell_id<-WhichCells(day6_8_atac)
cell_id_day6<-c(cell_id[grep("-1",cell_id)])
cell_id_day8<-c(cell_id[grep("-2",cell_id)])
cell_stat<-data.frame(cell=c(cell_id_day6,cell_id_day8),stat=c(rep("day6",length(cell_id_day6)),rep("day8",length(cell_id_day8))))
rownames(cell_stat)<-cell_stat$cell
day6_8_atac@meta.data$stat<-cell_stat[rownames(day6_8_atac@meta.data),]$stat
day6_atac<-day6_8_atac[,cell_id_day6]
day8_atac<-day6_8_atac[,cell_id_day8]
library(ggplot2)
plot1<-DimPlot(object = day6_atac, label = F,reduction = 'umap') +labs(title = "day6 scATAC")
plot2<-DimPlot(object = day8_atac, label = F,reduction = 'umap') +labs(title = "day8 scATAC")
CombinePlots(list(plot1,plot2),ncol = 1)

##figure4 c

###cicero
library(cicero)
DefaultAssay(day6_8_atac) <-'day6_8_peaks'
day6_8_atac.cds <- as.cell_data_set(x = day6_8_atac)
day6_8_atac.cicero <- make_cicero_cds(day6_8_atac.cds, reduced_coordinates = reducedDims(day6_8_atac.cds)$UMAP)
day6_8_atac_genome <- seqlengths(day6_8_atac)
# use chromosome 1 to save some time
# omit this step to run on the whole genome
day6_8_atac_genome <- day6_8_atac_genome[4]

# convert chromosome sizes to a dataframe
day6_8_atac_genome.df <- data.frame("chr" = names(day6_8_atac_genome), "length" = day6_8_atac_genome)

# run cicero
conns <- run_cicero(day6_8_atac.cicero, genomic_coords = day6_8_atac_genome.df)
ccans <- generate_ccans(conns)

links <- ConnectionsToLinks(conns = conns, ccans = ccans,threshold = 1)
Links(day6_8_atac) <- links
day6_8_gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(day6_8_gene.coords) <- 'UCSC'
day6_8_genebody.coords <- keepStandardChromosomes(day6_8_gene.coords, pruning.mode = 'coarse')
day6_8_genebodyandpromoter.coords <- Extend(x = day6_8_gene.coords, upstream = 2000, downstream = 0)

region2 <- GRangesToString(subset(day6_8_gene.coords, symbol=="Pax7"))
CoveragePlot(day6_8_atac, region = region2)
###separate condition
cell_id<-WhichCells(day6_8_atac)
cell_id_day6<-c(cell_id[grep("-1",cell_id)])
cell_id_day8<-c(cell_id[grep("-2",cell_id)])

day6_atac<-day6_8_atac[,cell_id_day6]
day8_atac<-day6_8_atac[,cell_id_day8]

#day6
library(SeuratWrappers)
DefaultAssay(day6_atac) <-'day6_8_peaks'
day6_atac.cds <- as.cell_data_set(x = day6_atac)

day6_atac.cicero <- make_cicero_cds(day6_atac.cds, reduced_coordinates = reducedDims(day6_atac.cds)$UMAP,k=500)
day6_atac_genome <- seqlengths(day6_atac)
# use chromosome 1 to save some time
# omit this step to run on the whole genome
day6_atac_genome_chr1 <- day6_atac_genome[1]
day6_atac_genome_chr4 <- day6_atac_genome[4]

# convert chromosome sizes to a dataframe
day6_atac_genome_chr1.df <- data.frame("chr" = names(day6_atac_genome_chr1), "length" = day6_atac_genome_chr1)
day6_atac_genome_chr4.df <- data.frame("chr" = names(day6_atac_genome_chr4), "length" = day6_atac_genome_chr4)

# run cicero
day6_atac_conns_chr1 <- run_cicero(day6_atac.cicero, genomic_coords = day6_atac_genome_chr1.df)
day6_atac_conns_chr1_sel <-subset(day6_atac_conns_chr1,coaccess>=0.3)
day6_atac_ccans_chr1 <- generate_ccans(day6_atac_conns_chr1_sel)

day6_atac_conns_chr4 <- run_cicero(day6_atac.cicero, genomic_coords = day6_atac_genome_chr4.df)


day6_atac_conns_chr4_sel <-subset(day6_atac_conns_chr4,coaccess>=0.3)
day6_atac_ccans_chr4 <- generate_ccans(day6_atac_conns_chr4_sel)

day6_atac_links_chr1 <- ConnectionsToLinks(conns = day6_atac_conns_chr1_sel, ccans = day6_atac_ccans_chr1)
day6_atac_links_chr4 <- ConnectionsToLinks(conns = day6_atac_conns_chr4_sel, ccans = day6_atac_ccans_chr4)

Links(day6_atac) <- c(day6_atac_links_chr1,day6_atac_links_chr4)
day6_8_gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(day6_8_gene.coords) <- 'UCSC'
day6_8_genebody.coords <- keepStandardChromosomes(day6_8_gene.coords, pruning.mode = 'coarse')
day6_8_genebodyandpromoter.coords <- Extend(x = day6_8_gene.coords, upstream = 2000, downstream = 0)

region1 <- GRangesToString(subset(day6_8_gene.coords, symbol=="Pax3"))
region2 <- GRangesToString(subset(day6_8_gene.coords, symbol=="Pax7"))
p1<-CoveragePlot(day6_atac, region = c(region1, region2),
                 extend.upstream = 2000,
                 extend.downstream = 2000,
                 idents = c("0","3","4"),
                 ymax = 19,
                 ncol = 2)
save(day6_atac_conns_chr1,day6_atac_conns_chr4,file="day6_cicero_chr1_chr4_v1_1_0.RData")

#day8
DefaultAssay(day8_atac) <-'day6_8_peaks'
day8_atac.cds <- as.cell_data_set(x = day8_atac)
day8_atac.cicero <- make_cicero_cds(day8_atac.cds, reduced_coordinates = reducedDims(day8_atac.cds)$UMAP,k=500)
day8_atac_genome <- seqlengths(day8_atac)
# use chromosome 1 to save some time
# omit this step to run on the whole genome
day8_atac_genome_chr1 <- day8_atac_genome[1]
day8_atac_genome_chr4 <- day8_atac_genome[4]

# convert chromosome sizes to a dataframe
day8_atac_genome_chr1.df <- data.frame("chr" = names(day8_atac_genome_chr1), "length" = day8_atac_genome_chr1)
day8_atac_genome_chr4.df <- data.frame("chr" = names(day8_atac_genome_chr4), "length" = day8_atac_genome_chr4)

# run cicero
day8_atac_conns_chr1 <- run_cicero(day8_atac.cicero, genomic_coords = day8_atac_genome_chr1.df)
day8_atac_conns_chr1_sel <-subset(day8_atac_conns_chr1,coaccess >= 0.7)
day8_atac_ccans_chr1 <- generate_ccans(day8_atac_conns_chr1_sel)

day8_atac_conns_chr4 <- run_cicero(day8_atac.cicero, genomic_coords = day8_atac_genome_chr4.df)
day8_atac_conns_chr4_sel <-subset(day8_atac_conns_chr4,coaccess >= 0.7)
day8_atac_ccans_chr4 <- generate_ccans(day8_atac_conns_chr4_sel)

day8_atac_links_chr1 <- ConnectionsToLinks(conns = day8_atac_conns_chr1_sel, ccans = day8_atac_ccans_chr1)
day8_atac_links_chr4 <- ConnectionsToLinks(conns = day8_atac_conns_chr4_sel, ccans = day8_atac_ccans_chr4)

Links(day8_atac) <- c(day8_atac_links_chr1,day8_atac_links_chr4)
day6_8_gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(day6_8_gene.coords) <- 'UCSC'
day6_8_genebody.coords <- keepStandardChromosomes(day6_8_gene.coords, pruning.mode = 'coarse')
day6_8_genebodyandpromoter.coords <- Extend(x = day6_8_gene.coords, upstream = 2000, downstream = 0)

region1 <- GRangesToString(subset(day6_8_gene.coords, symbol=="Pax3"))
region2 <- GRangesToString(subset(day6_8_gene.coords, symbol=="Pax7"))

p2<-CoveragePlot(day8_atac, region = c(region1, region2),
                 extend.upstream = 2000,
                 extend.downstream = 2000,
                 idents = c("1","3","4"),
                 ymax = 19,
                 ncol = 2)



CombinePlots(list(p1,p2),ncol = 1)
