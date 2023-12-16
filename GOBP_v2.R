library(escape)
library(ggplot2)
gene.sets1 <- getGeneSets(library = "C5", gene.sets = c('HP_RETINAL_NEOVASCULARIZATION',
            'GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS',
            'GOBP_BLOOD_VESSEL_MORPHOGENESIS','GOBP_RETINAL_BLOOD_VESSEL_MORPHOGENESIS',
            'GOBP_GLYCOLYTIC_PROCESS','GOBP_CELLULAR_LIPID_METABOLIC_PROCESS',
            'GOBP_LIPID_METABOLIC_PROCESS','GOBP_GLUCOSE_METABOLIC_PROCESS'), species = 'Mus musculus')
ES <- enrichIt(obj = sunye.combined,
gene.sets = gene.sets1,
groups = 3000)

sunye.combined <- AddMetaData(sunye.combined, ES)
sunye.combined <- SetIdent(sunye.combined, value = 'labels')

ES2 <- data.frame(sunye.combined[[]], Idents(sunye.combined))

colnames(ES2)[ncol(ES2)] <- "cluster"
head(ES2)
View(gene.sets1)

metadata_df <- data.frame(
  labels = sunye.combined$labels,
  groups = sunye.combined$group,
  GOBP = sunye.combined$HP_RETINAL_NEOVASCULARIZATION
)

average_expression <- aggregate(GOBP ~ labels + groups, data = metadata_df, FUN = mean)
# Displaying the table
print(average_expression)

average_expression$groups <- factor(average_expression$groups, levels = c('OIR_KO','OIR_WT','Normal'))

ggplot(average_expression, aes(x = labels, y = groups, fill= GOBP)) +
  geom_tile(color = "black",
            lwd = 1,linetype = 1) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "#FF0000", midpoint =-1150) + theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
