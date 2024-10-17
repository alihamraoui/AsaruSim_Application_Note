rename_mtx_with_gene_names <- function(matx.LR){
  
  ensemblsIDS <- rownames(matx.LR)
  genes_name <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensemblsIDS, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  
  known_genes <- length(rownames(matx.LR[which(rownames(matx.LR) %in% genes_name$GENEID),]))
  print(paste("nb of Unknown gene IDs:", length(rownames(matx.LR)) - known_genes))
  
  matx.LR <-matx.LR[which(rownames(matx.LR) %in% genes_name$GENEID),]
  matx.LR$GENEID <- rownames(matx.LR)
  
  matx.LR <- inner_join(matx.LR, genes_name, by = "GENEID")
  matx.LR <- matx.LR[,-which(colnames(matx.LR) == "GENEID")]
  
  dup_names <- matx.LR[duplicated(matx.LR$SYMBOL),'SYMBOL']
  if(length(dup_names)>0){
    dup_mtx <- matx.LR[which(matx.LR$SYMBOL %in% dup_names),]
    matx.LR <-  matx.LR[-which(matx.LR$SYMBOL %in% dup_names),]
    dup_mtx <- dup_mtx %>% aggregate(. ~ SYMBOL, FUN=sum)
    matx.LR <- rbind(matx.LR, dup_mtx)
  }
  rownames(matx.LR) <- matx.LR$SYMBOL
  matx.LR <- matx.LR[,-which(colnames(matx.LR) == "SYMBOL")]
  return(matx.LR)
}

clustuerer <- function (ScObject, min.features, nfeatures, dim) {
  ScObject = Seurat::NormalizeData(ScObject, normalization.method = "LogNormalize", verbose = F)
  ScObject = Seurat::FindVariableFeatures(ScObject, selection.method = "vst", nfeatures = nfeatures, verbose = F)
  all.genes <- rownames(ScObject)
  ScObject <- Seurat::ScaleData(ScObject, features = all.genes, verbose = F)
  ScObject <- Seurat::RunPCA(ScObject, features = Seurat::VariableFeatures(object = ScObject), verbose = F)
  ScObject <- Seurat::FindNeighbors(ScObject, dims = 1:dim, verbose = F)
  ScObject <- Seurat::FindClusters(ScObject, resolution = 0.1, verbose = F)
  ScObject <- Seurat::RunTSNE(ScObject, dims = 1:dim, verbose = F)
  ScObject <- Seurat::RunUMAP(ScObject, dims = 1:dim, verbose = F)
}

auto.annotation <- function(mtx, obj) {
  
  ref <- celldex::BlueprintEncodeData()
  
  pred <- SingleR(test = mtx, ref = ref, labels = ref$label.main)
  pred <- data.frame(BC = rownames(pred), cell_type = pred$labels)
  rownames(pred) <- pred$BC
  
  id_clusters <- data.frame("clusters"=obj$seurat_clusters)
  id_df <- cbind(id_clusters, pred[rownames(id_clusters), ])
  
  obj$cell_type <- id_df$cell_type
  
  majority_cell_type <- tapply(id_df$cell_type, id_df$clusters, function(x) {
    names(which.max(table(x)))
  })
  
  obj$corrected_cell_type <- majority_cell_type[id_df$clusters]
  
  return(obj)
}

integrate <- function(sim.Obj, real.Obj){
  obj <- merge(sim.Obj, real.Obj)
  obj <- Seurat::NormalizeData(obj, verbose = F)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  obj <- Seurat::ScaleData(obj, verbose = F)
  obj <- Seurat::RunPCA(obj, npcs =15, verbose = F)
  obj <- Seurat::IntegrateLayers(object = obj, method = Seurat::CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = F)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- Seurat::RunUMAP(obj, dims = 1:15, reduction = "integrated.cca", verbose = F)
  obj <- Seurat::RunTSNE(obj, dims = 1:15, reduction = "integrated.cca", verbose = F)
}

calc_meanLISI <- function(obj, reduction, ident="orig.ident"){
  if (ident == "cell_type"){
    meta_data <- data.frame(colnames(obj), ident=obj@meta.data$cell_type)
  }else{
    meta_data <- data.frame(colnames(obj), ident=obj@meta.data$orig.ident)
  }
  cord.umap <- Embeddings(obj, reduction = reduction)[, 1:2]
  lisi <- lisi::compute_lisi(cord.umap, meta_data, c("ident"))
  medianLISI <- median(lisi$ident)
}

DimPlot.integ <- function(Sobj, reduction) {
  lisi <- calc_meanLISI(Sobj, "umap")
  DimPlot(Sobj, reduction = reduction,  group.by = c("orig.ident"), cols = c('#FF7F0E', '#9467BD')) + 
    ggtitle(paste('miLISI  = ', round(lisi,2)-1)) + 
    theme(aspect.ratio = 1)}