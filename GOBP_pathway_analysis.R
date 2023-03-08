sunye_normal$state <- 'Normal'
sunye_OIRKO$state <- 'KOOIR'
sunye_OIRWT$state <- 'WTOIR'
sunye <- merge(sunye_normal, y = c(sunye_OIRWT, sunye_OIRKO))
library(escape)
gene.sets <- getGeneSets(library = "C5", species = 'Mus musculus', gene.sets = c('GOBP_PHAGOCYTOSIS','GOBP_MACROPHAGE_ACTIVATION','GOBP_MACROPHAGE_DIFFERENTIATION'))
ES <- enrichIt(obj = sunye,
gene.sets = gene.sets,
groups = 1000, cores = 4)
sunye <- AddMetaData(sunye, ES)
sunye <- SetIdent(sunye, value = 'labels')
ES2 <- data.frame(sunye [[]], Idents(sunye ))
colnames(ES2)[ncol(ES2)] <- "cluster"
View(gene.sets)
ridgeEnrichment(ES2, gene.set = 'GOBP_MACROPHAGE_ACTIVATION', group = "labels", add.rug = TRUE, facet = 'state')
p <- ES2 %>%
ggplot( aes(x=GOBP_PHAGOCYTOSIS, fill=state)) +
geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 30) +
scale_fill_manual(values=c("#69b3a2", "#404080")) +
theme_ipsum() +
labs(fill="")
