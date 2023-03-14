sunye <- readRDS("G://Sun_Ye/Processed/P13_0.5_3.RDS") #read the preprocessed file
sunye_normal <- subset(sunye, subset = group == 'Normal') #subset by conditions
head(sunye_normal)
DimPlot(sunye_normal, label = TRUE, label.box = TRUE)
sunye_normal$labels <- sunye_normal@active.ident
sunye_OIRKO <- subset(sunye, subset = group == 'OIR_KO')
sunye_OIRKO$labels <- sunye_OIRKO@active.ident
sunye_OIRWT <- subset(sunye, subset = group == 'OIR_WT')
sunye_OIRWT$labels <- sunye_OIRWT@active.ident #save labels in a separated column

library(CellChat) #perform CellChat analysis
cellchat <- createCellChat(object = sunye_normal, group.by = "labels") #grouped by cell populations
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 3) 
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_aggregate(cellchat, signaling = c('SPP1'), layout = "chord") #SPP1 as example, can dig into L-R pairs, contribution and specific L-R pair chord diagram
