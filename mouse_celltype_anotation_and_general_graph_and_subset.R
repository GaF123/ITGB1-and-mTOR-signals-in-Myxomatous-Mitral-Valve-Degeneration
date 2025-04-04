##marker genes expression
FeaturePlot(MV2, features = c("Col1a1","Dpt","Lum","Postn","Pdgfra","H2B-EGFP","Dcn","Prrx1"),reduction = "umap")

FeaturePlot(MV2, features = c("Pecam1", "Cdh5", "Vwf", "Icam2", "Kdr", "Sox18"),reduction = "umap")

FeaturePlot(MV2, features = c("Dct", "Pmel", "Tyrp1",  "Mlana","Mlph","Mitf"), reduction = "umap")

FeaturePlot(MV2, features = c("Csf1r","Cd68","Fcgr1","Itgam","Ccr2","Cd14","Cd115","Tnf","Nos2","Cd163","Arg1"),reduction = "umap")

FeaturePlot(MV2, features = c("Ptprc", "Cd52", "Cd48", "Cd53", "Itgb2", "Vav1"),reduction = "umap")

FeaturePlot(MV2, features = c("Mki67", "Ccnb2", "Top2a","Birc5","Cenpf", "Cdk1"),reduction =  "umap")

FeaturePlot(MV2, features = c("Cd3e", "Cd3g", "Cd3d","Cd4","Cd8a","Cd8b","Cd25","Foxp3","Tbx21","Ifng","Gata3","Il4","Rorc","Il17a"))

FeaturePlot(MV2, features = c("Cd19", "Cd79a", "Cd79b","Cd20","Sdc1","Prdm1","Cd138","Cd138","Prdm1"))

FeaturePlot(MV2, features = c("Klrb1b", "Ncr1", "Gzmb","Ncam1","Kir2d","Nkg2d"))

FeaturePlot(MV2, features = c("Cxcr2", "S100a9", "S100a8","Elane","Mpo","Cd16"))

FeaturePlot(MV2, features = c("Cd200r3", "Ms4a2","Kit","Fce1a"))

FeaturePlot(MV2, features = c("Myh11", "Acta2", "Itga8", "Myocd","Cnn1","Tagln","Smtn"))

FeaturePlot(MV2, features = c("Flt3","Itgax","Zbtb46","Lamp3","Cx3cr1","Itgam","Cd11c","Itgax","Cd1c"))

FeaturePlot(MV2, features = c("Hba1", "Hba2", "Hbb", "Hbg1", "Hbg2"))


##cluster definition:
levels(MV2)
new.cluster.ids <- c("Fibroblast","Fibroblast", "Fibroblast", "Leukocyte", "Endothelium","Fibroblast", "Leukocyte", "Leukocyte","Melanocyte","Leukocyte")
names(new.cluster.ids) <- levels(MV2)
new.cluster.ids
MV2 <- RenameIdents(MV2, new.cluster.ids)
DimPlot(MV2, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(MV2, reduction = "umap", label = FALSE, pt.size = 0.5) 

######## Add the new cluster names as a metadata column in your Seurat object
MV2[["cell_type"]] <- Idents(MV2)




##cell_type_percentage
cell.prop<-as.data.frame(prop.table(table(Idents(MV2), MV2$orig.ident)))
colnames(cell.prop)<-c("clusters","orig.ident","proportion")
pdf(file = "cell_type_percentage.pdf", width=6.5, height=6)
ggplot(cell.prop,aes(orig.ident,proportion,fill=clusters))+
  geom_bar(stat="identity",position="fill")+ggtitle("Proportion of cells in each culster")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'))+guides(fill=guide_legend(title=NULL))


## extract "Fibroblast",  "Endothelium", "Leukocyte", "Melanocyte"
Fibroblast <- subset(MV2, idents = c("Fibroblast"))
Fibroblast
saveRDS(Fibroblast,file = "Fibroblast.rds")


Leukocyte <- subset(MV2, idents = c("Leukocyte"))
Leukocyte 
saveRDS(Leukocyte,file = "Leukocyte.rds")


Endothelium <- subset(MV2, idents = "Endothelium")
Endothelium
saveRDS(Endothelium,file = "Endothelium.rds")


Melanocyte <- subset(MV2, idents = "Melanocyte")
Melanocyte
saveRDS(Melanocyte,file = "Melanocyte.rds")

