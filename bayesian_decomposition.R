library(Seurat)
library(SeuratDisk)

sunye_norm <- LoadH5Seurat('C://Users/rodri/Downloads/sunyennorm.h5Seurat', misc = F)
sunye_WTOIR <- LoadH5Seurat('C://Users/rodri/Downloads/sunyeOIRWT.h5Seurat', misc = F)
sunye_KOOIR <- LoadH5Seurat('C://Users/rodri/Downloads/sunyeOIRKO.h5Seurat', misc = F)

sunye_norm$labels <- sunye_norm@active.ident
sunye.combined <- merge(sunye_norm, y = sunye_WTOIR)
sunye.combined <- merge(sunye.combined, y = sunye_KOOIR)

library(scDC)
library(broom.mixed)
exprsMat <- GetAssayData(object = sunye.combined, assay = "RNA", slot = "data")
cellTypes <- sunye.combined$labels
subject <- sunye.combined$group
cond <- sunye.combined$group
dim(exprsMat)
table(subject, cellTypes)
table(cond, cellTypes)
res_scDC_noClust <- scDC_noClustering(cellTypes, subject, calCI = TRUE, 
                                      calCI_method = c("percentile", "BCa", "multinom"),
                                      nboot = 50)
barplotCI(res_scDC_noClust, c("Normal","WTOIR", 'KOOIR'))
densityCI(res_scDC_noClust, c("Normal", 'Normal',"Normal", 'Normal',"Normal", 'Normal',
                              "Normal", 'Normal',"Normal", 'Normal', 'Normal',
                              "WTOIR","WTOIR","WTOIR","WTOIR","WTOIR","WTOIR","WTOIR",
                              "WTOIR","WTOIR","WTOIR","WTOIR",'KOOIR','KOOIR','KOOIR',
                              'KOOIR','KOOIR','KOOIR','KOOIR','KOOIR','KOOIR','KOOIR','KOOIR'))
res_GLM <- fitGLM(res_scDC_noClust, c("Normal", 'Normal',"Normal", 'Normal',"Normal", 'Normal',
                                      "Normal", 'Normal',"Normal", 'Normal', 'Normal',
                                      "WTOIR","WTOIR","WTOIR","WTOIR","WTOIR","WTOIR","WTOIR",
                                      "WTOIR","WTOIR","WTOIR","WTOIR",'KOOIR','KOOIR','KOOIR',
                                      'KOOIR','KOOIR','KOOIR','KOOIR','KOOIR','KOOIR','KOOIR','KOOIR'), 
                  pairwise = T)
summary(res_GLM$pool_res_fixed)
summary(res_GLM$pool_res_random)
