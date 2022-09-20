library(CellChat)
library(ggplot2)
library(ggalluvial)
options(stringsAsFactors = FALSE)


Load GBM R-object
Idents(object= GBM) <- GBM@meta.data$cellType
data.input <- GetAssayData(GBM, assay = "RNA", slot = "data")
labels <- Idents(GBM)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(data = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the objec
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
options(future.globals.maxSize= 1691289600)
future::plan("multiprocess", workers = 6) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
netVisual_circle(cellchat@net$count, vertex.size = groupSize, weight.scale = T, label.edge= F, edge.label.cex = 0.8, vertex.label.cex = 1)
vertex.receiver = seq(1,3)
##output all the significant pathway
cellchat@netP$pathways
pathways.show <- "TGFb"
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
netVisual_aggregate(cellchat, signaling = c("TGFb"), layout = "circle", vertex.size = groupSize)
netAnalysis_contribution(cellchat, signaling = pathways.show)
netVisual_aggregate(cellchat, signaling = c("IL6"), layout = "circle", vertex.size = groupSize)
netVisual_aggregate(cellchat, signaling = c("CSF"), layout = "circle", vertex.size = groupSize)