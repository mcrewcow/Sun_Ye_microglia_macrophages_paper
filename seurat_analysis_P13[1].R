#Guoshuai Cai's pipeline for scRNA-seq data processing, analysis and visualization for results shown in Figs. 3C-E, 4A,B and 5A,B.

library(dplyr)
library(Seurat)
library(cowplot)

batch<-"P13"
res<-0.5

set.seed(123)

dir<-"."

dir.in<-"/Users/gcai/Dropbox/work/Collabration/YeSun/Single cell seq Data/Project 1"

dir.in<-file.path(dir.in,"P13N and P13H")

sample.list<-list(Normal="Normal WT Pool84-4_3",
	OIR_WT="OIR WT Pool84-4_1",
	OIR_KO="OIR KO Pool84-4_2")

filename<-batch

dir.out<-file.path(dir,"Seurat")

if (file.exists(file.path(dir.out))){
} else {
    	dir.create(file.path(dir.out))
}

dir.out<-file.path(dir.out,filename)

if (file.exists(file.path(dir.out))){
} else {
    	dir.create(file.path(dir.out))
}

sample.vec<-names(sample.list)

obj.list<-list()

for(sample.i in sample.vec){

	dir.i<-file.path(dir.in,sample.list[[sample.i]])

	data.i <- Read10X(data.dir = file.path(dir.i,"filtered_feature_bc_matrix"))
	obj.i <- CreateSeuratObject(counts = data.i, project = sample.i, min.cells = 3, min.features = 200)
	obj.i$group <- sample.i

	obj.list[[sample.i]]<-obj.i
}

obj <- Reduce(merge, obj.list)

#normalization
obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)

#FindVariableFeatures
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

#store mitochondrial percentage in object meta data
obj <- PercentageFeatureSet(obj, pattern = "^MT-|^Mt-|^mt-", col.name = "percent.mt")


obj@active.assay <- "RNA"

pdf(file = file.path(dir.out,paste("FeatureScatter_",filename,".pdf",sep="")))
	par(mfrow = c(1, 2))
	FeatureScatter(object = obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
	FeatureScatter(object = obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

#filtering
if(batch=="P13"){
	obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 8)
}

# Run the standard workflow for visualization and clustering
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

#Determine the ‘dimensionality’ of the dataset
obj <- JackStraw(obj, num.replicate = 100)
obj <- ScoreJackStraw(obj, dims = 1:20)

pdf(file = file.path(dir.out,paste("JackStrawPlot_",filename,".pdf",sep="")))
	JackStrawPlot(obj, dims = 1:15)
dev.off()

####dimensionality determination strategy 1
pdf(file = file.path(dir.out,paste("ElbowPlot_",filename,".pdf",sep="")))
	ElbowPlot(obj)
dev.off()

# t-SNE and Clustering
obj <- RunTSNE(obj, reduction = "pca", dims.use = 1:20)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:20)
obj <- FindClusters(obj, resolution = res)

save(obj,file=file.path(dir.out,paste0(filename,"_",res,".obj.RData")))

#assign cell types
obj$orig.ident<-Idents(obj)
Idents(obj)<-gsub("_.+","",Idents(obj))

if(batch=="P13"){
	obj <- RenameIdents(obj,   `0` = "microglia_1", 
    `1` = "microglia_1", `2` = "macrophage_1",`3` = "muller",`4` = "macrophage_2", `5` = "microglia_2", `6` = "B_1", `7` = "microglia_3",`8` = "muller", `9` = "monocyte",
    `10` = "microglia_3",`11` = "B_2",`12`="microglia_4",`13`="T",`14`="unknown_photoreceptor",`15` = "granulocyte",`16` = "stromal",`17`="endothelial_1",`18`="endothelial_2",`19` = "endothelial_2")
}

#filter CD45- cells and reclustering
obj<-subset(obj, idents=c("endothelial_1","endothelial_2","stromal","unknown_photoreceptor","muller"),invert=TRUE)

obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunTSNE(obj, reduction = "pca", dims.use = 1:20)

#DE between groups, for Figure 5B
cell.type<-unique(Idents(obj))

Idents(obj)<-paste(Idents(obj), obj$group, sep = ".")

group<-unique(obj$group)

for(cell.i in cell.type){
	for(m in length(group):2){
		for(n in (m-1):1){
			id<-paste(cell.i,group[c(m,n)],sep=".")
			DE<- try(FindMarkers(obj, ident.1 = id[1], ident.2 = id[2], verbose = FALSE, min.pct=0,logfc.threshold=0),silent=TRUE)
			if(class(DE)=="try-error"){
				DE<-DE[1]
			}else{
				write.csv(DE, file = file.path(dir.out,paste0(cell.i,".",paste(group[c(m,n)],collapse="_vs_"),".csv")),quote=FALSE)
			}
		}
	}
}

saveRDS(obj,file=file.path(dir.out,paste0(filename,"_",res,".RDS")))

#tSNE plot, for Fig. 2C
p1 <- DimPlot(obj, reduction = "tsne", group.by = "group", label = FALSE, pt.size=1.5)
p2 <- DimPlot(obj, reduction = "tsne", label = FALSE, pt.size=1)
pdf(file = file.path(dir.out,paste("tsne_",filename,"_",res,".pdf",sep="")),width=6.5, height=5)
	plot(p1)
	plot(p2)
dev.off()

#marker dot plot, for Fig. 2D
marker.vec<-c("P2ry12","Olfml3","Cx3cr1","Pf4","Ccl7","Ms4a7","Ms4a4a","Lyve1","Clec4n","Cyr61",
	"Rlbp1","Clu","Crym","Car2","Lgals3","Anxa2","Lpl","Hilpda","Cdkn1a","Cd79a","Ly6d","Ms4a1","Iglc2","Cd79b",
	"Birc5","Pclaf","Top2a","Ube2c","Ifitm6","Pglyrp1","Cybb","Ifitm3","Ifitm2","Irf7","Trim30a","Isg15",
	"Ms4a4b","Cd3d","Trbc2","Cd3e","Cd3g","S100a9","S100a8","Stfa1","Gm5483","Stfa2l1")
pdf(file = file.path(dir.out.plot,paste("DotPlot_celltype_marker_",filename,".pdf",sep="")),width=16,height=5)
	DotPlot(obj, features = marker.vec, dot.scale = 6) + RotatedAxis()
dev.off()

#tSNE plot,  for Fig. 2E
p1 <- DimPlot(obj, reduction = "tsne", split.by="group",label = FALSE, pt.size = 1.5)
pdf(file = file.path(dir.out,paste("tsne_split_",filename,"_",res,".pdf",sep="")),width=15, height=5)
	plot(p1)
dev.off()

#highlighted cell type plot, for Fig. 4B
require(scales)
cols<-hue_pal()(length(levels(Idents(obj))))
pdf(file = file.path(dir.out,paste("tsne_split_highlight_",filename,"_",res,".pdf",sep="")),width=15, height=5)
	for(type.i in c("microglia_1","microglia_2","microglia_3","macrophage_2")){
		cells.i <- WhichCells(obj, idents = type.i)
		col.i<-cols[which(levels(Idents(obj))==type.i)]
		p1 <- DimPlot(obj, reduction = "tsne", split.by="group",label = FALSE, pt.size = 1.5, cells.highlight=cells.i, cols.highlight=col.i,sizes.highlight = 2)
			plot(p1)
	}
dev.off()

#Spp1 feature plot, for Fig. 5B
pdf(file = file.path(dir.out,paste("FeaturePlot_target_",filename,".pdf",sep="")),width=10,height=10)
	FeaturePlot(obj, features = c("Spp1"), split.by = "group", reduction = "tsne", cols = c("grey", "red"), pt.size=1, order=TRUE)
dev.off()
