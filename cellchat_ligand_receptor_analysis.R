library(CellChat)
library(patchwork)
library(tidyverse)
library(Seurat)
library(ggplot2)
options(stringsAsFactors = FALSE)

cellchat <- createCellChat(object = mgR, meta = mgR@meta.data, group.by = "celltype_all")
cellchat

CellChatDB <- CellChatDB.mouse  
showDatabaseCategory(CellChatDB)


dplyr::glimpse(CellChatDB$interaction)



# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat@DB <- CellChatDB



cellchat <- subsetData(cellchat) 
##future::plan("multiprocess", workers = 4)
future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)



cellchat <- computeCommunProb(cellchat) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)




groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot1.pdf",width =6, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()







mat <- cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot2.pdf",width =6, height=6)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()





pathways.show <- c("CXCL")  

#levels(cellchat@idents)   
#[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
vertex.receiver = c(1,2,3,4,5) 
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot3.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.label.cex = 0.5,
                    vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot4.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

pathways.show <- c("TGFB")  
#levels(cellchat@idents)   
#[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
vertex.receiver = c(1,2,3,4,5) 
par(mfrow = c(1,2), xpd=TRUE)
pdf("plotTGFB.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.label.cex = 0.5,
                    vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot41.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

pathways.show <- c("WNT")
#levels(cellchat@idents)   
#[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
vertex.receiver = c(1,2,3,4,5) 
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot5.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.label.cex = 0.5,
                    vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()



par(mfrow = c(1,2), xpd=TRUE)
pdf("plot6.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
dev.off()

pathways.show <- c("CXCL")  
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot8.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()


# Heatmap
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot9.pdf",width =10, height=10)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
#> Do heatmap based on a single object

par(mfrow = c(1,2), xpd=TRUE)
pdf("plot10.pdf",width =10, height=10)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()



pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot11.pdf",width =10, height=10)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()
#> [[1]]
# Circle plot
pdf("plot12.pdf",width =10, height=10)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot13.pdf",width =10, height=10)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()




# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 5, units = 'in', dpi = 300)
}




# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot15.pdf",width =6, height=50)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
dev.off()
#> Comparing communications on a single object




# show all the significant interactions (L-R pairs) associated with certain signaling pathways
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot16.pdf",width =6, height=5)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
dev.off()
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot17.pdf",width =6, height=5)
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
dev.off()
#> Comparing communications on a single object



# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot18.pdf",width =20, height=20)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 1, legend.pos.y = 30)
dev.off()


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot19.pdf",width =20, height=20)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
dev.off()



# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot20.pdf",width =20, height=20)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf("plot21.pdf",width =20, height=20)
plotGeneExpression(cellchat, signaling = "SPP1")
dev.off()
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf("plot22.pdf",width =20, height=20)
plotGeneExpression(cellchat, signaling = "CXCL")
dev.off()
dev.off()


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot23.pdf",width =20, height=20)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()
dev.off()



# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot24.pdf",width =10, height=10)
gg1 + gg2
dev.off()
dev.off()


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
par(mfrow = c(1,2), xpd=TRUE,label.size = 0.001)
par(mar = c(1,1,1,1)+0.5)
pdf("plot25.pdf",width =60, height=60)
ht1 + ht2
dev.off()
dev.off()



# Signaling role analysis on the cell-cell communication networks of interest
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot26.pdf",width =10, height=10)
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))
dev.off()



# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))


selectK(cellchat, pattern = "outgoing")


nPatterns = 3
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot27.pdf",width =10, height=10)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()


# river plot
library(ggalluvial)
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot28.pdf",width =10, height=10)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function


# dot plot
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot29.pdf",width =10, height=10)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()
dev.off()



selectK(cellchat, pattern = "incoming")


nPatterns = 4
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot30.pdf",width =10, height=10)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()
dev.off()


# river plot
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot31.pdf",width =10, height=10)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function


# dot plot
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot32.pdf",width =10, height=10)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()
dev.off()



cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot33.pdf",width =10, height=10)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()
dev.off()



cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)     


netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


saveRDS(cellchat, file = "cellchat_mouse_mgR.rds")







library(CellChat)
library(patchwork)
library(tidyverse)
library(Seurat)
library(ggplot2)
options(stringsAsFactors = FALSE)


head(WT@meta.data)

cellchat <- createCellChat(object = WT, meta = WT@meta.data, group.by = "celltype_all")
cellchat

CellChatDB <- CellChatDB.mouse  
showDatabaseCategory(CellChatDB)


dplyr::glimpse(CellChatDB$interaction)



# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat@DB <- CellChatDB



cellchat <- subsetData(cellchat) 
##future::plan("multiprocess", workers = 4)
future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)



setwd("/Users/eaiscui/Desktop/cellchat/results/wt/")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot1.pdf",width =6, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()







mat <- cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot2.pdf",width =6, height=6)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()





pathways.show <- c("CXCL")  

#levels(cellchat@idents)   
#[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
vertex.receiver = c(1,2,3,4,5) 
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot3.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.label.cex = 0.5,
                    vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot4.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()



pathways.show <- c("WNT")
#levels(cellchat@idents)   
#[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
vertex.receiver = c(1,2,3,4,5) 
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot5.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.label.cex = 0.5,
                    vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf("plot6.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
dev.off()

pathways.show <- c("CXCL")  
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot8.pdf",width =10, height=10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()


# Heatmap
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot9.pdf",width =10, height=10)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
#> Do heatmap based on a single object

par(mfrow = c(1,2), xpd=TRUE)
pdf("plot10.pdf",width =10, height=10)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
dev.off()


pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot11.pdf",width =10, height=10)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()
#> [[1]]
# Circle plot
pdf("plot12.pdf",width =10, height=10)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot13.pdf",width =10, height=10)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()




# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 5, units = 'in', dpi = 300)
}




# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot15.pdf",width =6, height=50)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
dev.off()
#> Comparing communications on a single object




# show all the significant interactions (L-R pairs) associated with certain signaling pathways
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot16.pdf",width =6, height=5)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
dev.off()
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot17.pdf",width =6, height=5)
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
dev.off()
#> Comparing communications on a single object



# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot18.pdf",width =20, height=20)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 1, legend.pos.y = 30)
dev.off()


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot19.pdf",width =20, height=20)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
dev.off()



# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot20.pdf",width =20, height=20)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf("plot21.pdf",width =20, height=20)
plotGeneExpression(cellchat, signaling = "SPP1")
dev.off()
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf("plot22.pdf",width =20, height=20)
plotGeneExpression(cellchat, signaling = "CXCL")
dev.off()
dev.off()


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot23.pdf",width =20, height=20)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()
dev.off()


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot24.pdf",width =10, height=10)
gg1 + gg2
dev.off()
dev.off()


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
par(mfrow = c(1,2), xpd=TRUE,label.size = 0.001)
par(mar = c(1,1,1,1)+0.5)
pdf("plot25.pdf",width =30, height=30)
ht1 + ht2
dev.off()
dev.off()



# Signaling role analysis on the cell-cell communication networks of interest
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot26.pdf",width =10, height=10)
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))
dev.off()


# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))


selectK(cellchat, pattern = "outgoing")


nPatterns = 3
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot27.pdf",width =10, height=10)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()


# river plot
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot28.pdf",width =10, height=10)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function


# dot plot
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot29.pdf",width =10, height=10)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()
dev.off()



selectK(cellchat, pattern = "incoming")


nPatterns = 4
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot30.pdf",width =10, height=10)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()
dev.off()


# river plot
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot31.pdf",width =10, height=10)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function


# dot plot
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot32.pdf",width =10, height=10)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()
dev.off()



cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
par(mfrow = c(1,2), xpd=TRUE)
pdf("plot33.pdf",width =10, height=10)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()
dev.off()



cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)     


netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


saveRDS(cellchat, file = "cellchat_mouse_WT.rds")



library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)

data.dir <- '/Users/eaiscui/Desktop/cellchat/results/merged'
dir.create(data.dir)
setwd(data.dir)


cellchat.WT <- readRDS("/Users/eaiscui/Desktop/cellchatRDS/cellchat_mouse_WT.rds")
cellchat.mgR <- readRDS("/Users/eaiscui/Desktop/cellchatRDS/cellchat_mouse_mgR.rds")

object.list <- list(WT = cellchat.WT, mgR = cellchat.mgR)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 
#>  555 signaling genes.
#>  7563 cells. 
#> CellChat analysis of single cell RNA-seq data!

par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchat1.pdf", width = 20, height =16)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()
dev.off()




pdf(file ="cellchat1.pdf", width = 20, height =16)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) +
  theme(axis.text = element_text(size = 40), # change the size as needed
        axis.title = element_text(size = 40)) # change the size as needed
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") +
  theme(axis.text = element_text(size = 40), # change the size as needed
        axis.title = element_text(size = 40)) # change the size as needed
gg1 + gg2
dev.off()
dev.off()










par(mfrow = c(1,2), xpd=TRUE, cex=2)
pdf(file ="cellchat2.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE, cex=2)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchat3.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()
dev.off()













par(mfrow = c(1,2), xpd=TRUE, cex=2)
pdf(file ="cellchat4.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE, cex=2)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE, cex=2)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
dev.off()



par(mfrow = c(1,2), xpd=TRUE, cex=1)
pdf(file ="cellchat5.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE, cex=1)
par(mfrow = c(1,2), xpd=TRUE, cex=1)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T)
dev.off()
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchat6.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE)
par(mfrow = c(1,2), xpd=TRUE)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", label.edge = T)
dev.off()
dev.off()
par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchat7.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T)
dev.off()
dev.off()

par(mfrow = c(1,2), xpd=TRUE, cex=2)
pdf(file ="cellchat8.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE, cex=2)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.off()
dev.off()

#####################
par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchat9.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "4-Fibroblast")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
dev.off()
dev.off()
dev.off()
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchat10.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "4-Fibroblast")
gg2 <- gg2 + coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(-0.05, 0.05)) # 根据您的数据设置适当的范围

#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()
dev.off()
######################


cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2




par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchat11.pdf", width = 20, height =16)
par(mfrow = c(1,2), xpd=TRUE)
rankSimilarity(cellchat, type = "functional")
dev.off()
dev.off()
#> Compute the distance of signaling networks between datasets 1 2



gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2



library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.10.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd", font.size = 2.5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd", font.size = 2.5)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:14),  comparison = c(1, 2), angle.x = 45, font.size = 2.5)
#> Comparing communications on a merged object







par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchatincreseadecres1.pdf", width = 20, height =16)
gg1 <- netVisual_bubble(cellchat, sources.use = c(11,5,3,6), targets.use = c(4),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T, font.size = 14)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(11,5,3,6), targets.use = c(4),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T, font.size = 9)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchatincreseadecres2.pdf", width = 20, height =16)
gg1 <- netVisual_bubble(cellchat, sources.use = c(11,5,3,6,4,7), targets.use = c(1),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T, font.size = 14)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(11,5,3,6,4,7), targets.use = c(1),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T, font.size = 10)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()
dev.off()



signaling.LSIncreased = gg1$data.



# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "C1041G"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "C1041G",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(11,5,3,6), targets.use = c(4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(11,5,3,6), targets.use = c(4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(11,5,3,6,4,7), targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(11,5,3,6,4,7), targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2




pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(11,5,3,6), targets.use = c(4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(11,5,3,6), targets.use = c(4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2



par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cellchatincreseadecres3.pdf", width = 20, height =16)
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(11,5,3,6,4,7), targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), font.size = 18)
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(11,5,3,6,4,7), targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), font.size = 14)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()
dev.off()



pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(11,5,3,6), targets.use = c(4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(11,5,3,6), targets.use = c(4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2



gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(11,5,3,6,4,7), targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(11,5,3,6,4,7), targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2







# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cord1.pdf", width = 50, height =50)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:14), targets.use = c(4), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(1:14), targets.use = c(4), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
dev.off()

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cord2.pdf", width = 50, height =50)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:14), targets.use = c(4), slot.name = 'net', net = net.up, lab.cex = 2.5, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(1:14), targets.use = c(4), slot.name = 'net', net = net.down, lab.cex = 2.5, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cord3.pdf", width = 70, height =50)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(11,5,3,6), targets.use = c(4), slot.name = 'net', net = net.up, lab.cex = 2, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(11,5,3,6), targets.use = c(4), slot.name = 'net', net = net.down, lab.cex = 4, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cord4.pdf", width = 70, height =50)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(11,5,3,6), targets.use = c(4), slot.name = 'net', net = net.up, lab.cex = 5, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(11,5,3,6), targets.use = c(4), slot.name = 'net', net = net.down, lab.cex = 5, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cord5.pdf", width = 70, height =50)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(11,5,3,6,4,7), targets.use = c(1), slot.name = 'net', net = net.up, lab.cex = 5, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(11,5,3,6,4,7), targets.use = c(1), slot.name = 'net', net = net.down, lab.cex = 5, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
dev.off()



# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'mouse')


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "C1041G")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
#> The default behaviour of split.by has changed.
#> Separate violin plots are now plotted side-by-side.
#> To restore the old behaviour of a single split violin,
#> set split.plot = TRUE.
#>       
#> This message will be shown once per session.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.



cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "C1041G")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)



cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "C1041G")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.10.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu", font.size = 3)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu", font.size = 3)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))






# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cord6.pdf", width = 20, height =16)
netVisual_chord_gene(object.list[[2]], sources.use = c(11,5,3,6), targets.use = c(4), slot.name = 'net', net = net.up, lab.cex = 2, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(11,5,3,6), targets.use = c(4), slot.name = 'net', net = net.down, lab.cex = 2, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
dev.off()



par(mfrow = c(1,2), xpd=TRUE)
pdf(file ="cord7.pdf", width = 20, height =16)
netVisual_chord_gene(object.list[[2]], sources.use = c(11,5,3,6,4,7), targets.use = c(1), slot.name = 'net', net = net.up, lab.cex = 2, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(11,5,3,6,4,7), targets.use = c(1), slot.name = 'net', net = net.down, lab.cex = 2, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
dev.off()


saveRDS(cellchat,file = "cellchatfinal.rds")

