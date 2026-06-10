library(data.table)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(SingleR)
library(RColorBrewer)
library(Matrix)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyr)


### STEP 1: Download and inspect files ###
raw_data_P1A <- Read10X(data.dir = "GSE161277_RAW/P1/adenoma")
raw_data_P1C <- Read10X(data.dir = "GSE161277_RAW/P1/carcinoma")
raw_data_P1N <- Read10X(data.dir = "GSE161277_RAW/P1/normal")

raw_data_P2A <- Read10X(data.dir = "GSE161277_RAW/P2/adenoma")
raw_data_P2C <- Read10X(data.dir = "GSE161277_RAW/P2/carcinoma")
raw_data_P2N <- Read10X(data.dir = "GSE161277_RAW/P2/normal")

# Create Seurat object:
seurat_obj_P1A <- CreateSeuratObject(counts = raw_data_P1A, project = "GSE161277", 
                                     min.cells = 3, min.features = 200)

seurat_obj_P1C <- CreateSeuratObject(counts = raw_data_P1C, project = "GSE161277",
                                     min.cells = 3, min.features = 200)

seurat_obj_P1N <- CreateSeuratObject(counts = raw_data_P1N, project = "GSE161277", 
                                     min.cells = 3,min.features = 200)


seurat_obj_P2A <- CreateSeuratObject(counts = raw_data_P2A, project = "GSE161277", 
                                     min.cells = 3, min.features = 200)

seurat_obj_P2C <- CreateSeuratObject(counts = raw_data_P2C, project = "GSE161277", 
                                     min.cells = 3, min.features = 200)

seurat_obj_P2N <- CreateSeuratObject(counts = raw_data_P2N, project = "GSE161277", 
                                     min.cells = 3, min.features = 200)

#metadata
seurat_obj_P1A$Patient <- "P1"
seurat_obj_P1A$condition <- "Adenoma"

seurat_obj_P1C$Patient <- "P1"
seurat_obj_P1C$condition <- "Carcinoma"

seurat_obj_P1N$Patient <- "P1"
seurat_obj_P1N$condition <- "Normal"


seurat_obj_P2A$Patient <- "P2"
seurat_obj_P2A$condition <- "Adenoma"

seurat_obj_P2C$Patient <- "P2"
seurat_obj_P2C$condition <- "Carcinoma"

seurat_obj_P2N$Patient <- "P2"
seurat_obj_P2N$condition <- "Normal"

#Merge
P1 <- merge(seurat_obj_P1N, y = list(seurat_obj_P1A, seurat_obj_P1C),
  add.cell.ids = c("P1_Normal", "P1_Adenoma", "P1_Carcinoma"))


P2 <- merge(seurat_obj_P2N, y = list(seurat_obj_P2A, seurat_obj_P2C),
  add.cell.ids = c("P2_Normal", "P2_Adenoma", "P2_Carcinoma"))


### STEP 2: Quality Control ###
#P1
P1[["percent.mt"]] <- PercentageFeatureSet(P1, pattern = "^MT-")

VP_P1 <- VlnPlot(P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               pt.size = 0.1, ncol = 3)

summary(P1$nFeature_RNA)
summary(P1$nCount_RNA)

FS_P1 <- FeatureScatter(P1, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position = "none") + ggtitle("Correlation = -0.17")

P1 <- subset(P1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 
             & percent.mt < 20)
             
#P2
P2[["percent.mt"]] <- PercentageFeatureSet(P2, pattern = "^MT-")


VP_P2 <- VlnPlot(P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, 
        ncol = 3)

summary(P2$nFeature_RNA)
summary(P2$nCount_RNA)

FS_P2 <- FeatureScatter(P2, feature1 = "nCount_RNA", feature2 = "percent.mt") +  
  theme(legend.position = "none") + ggtitle("Correlation = -0.07")


P2 <- subset(P2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 &
               percent.mt < 20)

### STEP 3: Normalization and Feature Selection ###
#Merge patient 1 and patient 2 in one list
patients <- list(P1 = P1, P2 = P2)

patients <- lapply(patients, function(obj) {obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)
obj})

P1 <- patients$P1
P2 <- patients$P2

var_plot_1 <- VariableFeaturePlot(P1)
LP_P1 <- LabelPoints(plot = var_plot_1, points = head(VariableFeatures(P1), 20), 
            repel = TRUE) + ggtitle("Patient 1")

var_plot_2 <- VariableFeaturePlot(P2)
LP_P2 <- LabelPoints(plot = var_plot_2, points = head(VariableFeatures(P2), 20), 
                    repel = TRUE) + ggtitle("Patient 2")

#Dimensionality Reduction & Clustering
ElbowPlot(P1, ndims = 50)
ElbowPlot(P2, ndims = 50)

# Number of PCs selected based on ElbowPlot inspection
P1 <- FindNeighbors(P1, dims = 1:20)
P2 <- FindNeighbors(P2, dims = 1:22)

P1 <- FindClusters(P1, resolution = 0.6)
P2 <- FindClusters(P2, resolution = 0.6)

# UMAP
P1 <- RunUMAP(P1, dims = 1:20)
UMAP_P1 <- DimPlot(P1, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP: Clustered Cells")

P2 <- RunUMAP(P2, dims = 1:22)
UMAP_P2 <- DimPlot(P2, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP: Clustered Cells")

#UMAP based on condition 
UMAP_P1con <- DimPlot(P1, group.by = "condition")
UMAP_P2con <- DimPlot(P2, group.by = "condition")


### STEP 4: Cell Type Annotation ###
#P1
DefaultAssay(P1) <- "RNA"
P1 <- JoinLayers(P1)

Idents(P1) <- "seurat_clusters"

markers_P1 <- FindAllMarkers(P1, assay = "RNA", only.pos = TRUE, min.pct = 0.1, 
                          logfc.threshold = 0.1)

top_markers_P1 <- markers_P1 %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

canonical_markers <- c(
  "EPCAM", "KRT8", "KRT18", "KRT19",   # epithelial
  "CD3D", "CD3E", "TRAC",              # T cells
  "MS4A1", "CD79A", "CD74",            # B cells
  "NKG7", "GNLY",                      # NK cells
  "LYZ", "S100A8", "S100A9", "LST1",   # myeloid / monocytes
  "COL1A1", "COL1A2", "DCN", "ACTA2", # fibroblasts / CAF
  "VWF", "PECAM1",                     # endothelial
  "MKI67", "TOP2A",                    # proliferating
  "MUC2", "TFF3", "KRT20"              # intestinal epithelial / goblet-like
)

TM_P1 <- DotPlot(P1, features = canonical_markers, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("Markers in clusters, Patient 1")

FeaturePlot(P1, features = c("EPCAM", "CD3D", "MS4A1", "LYZ", "COL1A1", "PECAM1"), 
            ncol = 3)

new_ids_P1 <- c(
  "0" = "T cells",
  "1" = "epithelial + intestinal",
  "2" = "epithelial",
  "3" = "epithelial",
  "4" = "T cells",
  "5" = "T cells",
  "6" = "Proliferating epithelial",
  "7" = "immune-like",
  "8" = "epithelial",
  "9" = "B cells",
  "10" = "epithelial",
  "11" = "epithelial",
  "12" = "B cells",
  "13" = "NK cells",
  "14" = "epithelial",
  "15" = "myeloid",
  "16" = "Goblet-like",
  "17" = "epithelial",
  "18" = "Goblet-like",
  "19" = "epithelial"
)

#Verifying cluster 7
markers_P1 %>%
  filter(cluster == 7) %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

table(P1$seurat_clusters)

FeaturePlot(P1, features = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"),
            ncol = 3)
VlnPlot(P1, features = c("CD3D", "TRAC", "NKG7", "GNLY"), group.by = "seurat_clusters")
VlnPlot(P1, features = c("CD4", "CD8A", "IL7R", "CCR7"), group.by = "seurat_clusters",
  pt.size = 0)

#Cluster Visualization by cell type
P1$manual_celltype <- plyr::mapvalues(x = as.character(P1$seurat_clusters),
  from = names(new_ids_P1), to = new_ids_P1)
Idents(P1) <- P1$manual_celltype

DP_P1 <- DimPlot(P1, group.by = "manual_celltype", label = TRUE, repel = TRUE) +
  ggtitle("Annotated cell types, Patient 1")


#P2
DefaultAssay(P2) <- "RNA"
P2 <- JoinLayers(P2)

Idents(P2) <- "seurat_clusters"

markers_P2 <- FindAllMarkers(P2, assay = "RNA", only.pos = TRUE, min.pct = 0.1, 
                          logfc.threshold = 0.1)

top_markers_P2 <- markers_P2 %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

TM_P2 <- DotPlot(P2, features = canonical_markers, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("Markers in clusters, Patient 2")

FeaturePlot(P2, features = c("EPCAM", "CD3D", "MS4A1", "LYZ", "COL1A1", "PECAM1"), 
            ncol = 3)

new_ids_P2 <- c(
  "0" = "T cells",
  "1" = "B cells",
  "2" = "myeloid",
  "3" = "epithelial",
  "4" = "T cells",
  "5" = "T cells",
  "6" = "epithelial",
  "7" = "Proliferating epithelial",
  "8" = "intestinal",
  "9" = "NK cells",
  "10" = "Proliferating epithelial",
  "11" = "intestinal",
  "12" = "intestinal",
  "13" = "B cells",
  "14" = "epithelial",
  "15" = "fibroblasts / CAF",
  "16" = "intestinal",
  "17" = "epithelial",
  "18" = "endothelial",
  "19" = "B cells"
)

#Cluster Visualization by cell type
P2$manual_celltype <- plyr::mapvalues(x = as.character(P2$seurat_clusters),
                                      from = names(new_ids_P2), to = new_ids_P2)
Idents(P2) <- P2$manual_celltype

DP_P2 <- DimPlot(P2, group.by = "manual_celltype", label = TRUE, repel = TRUE) +
  ggtitle("Annotated cell types, patient 2")

table(P1$manual_celltype, P1$condition)
table(P2$manual_celltype, P2$condition)

#Celltypes differences based on conditions
DP_P1con <- DimPlot(P1, group.by = "manual_celltype", split.by = "condition")
DP_P2con <- DimPlot(P2, group.by = "manual_celltype", split.by = "condition")

run_celltype_diff <- function(seurat_obj, patient_label) {
  meta <- as.data.frame(seurat_obj@meta.data)
  celltype_counts <- meta %>% 
    dplyr::count(condition, manual_celltype) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(percent = n / sum(n) * 100)
  plot <- ggplot(celltype_counts, aes(x = condition, y = percent, fill = manual_celltype)) +
    geom_col(position = "fill") + scale_fill_brewer(palette = "Set3") + theme_minimal() + 
    labs(title = paste("Cell type composition by condition,", patient_label), 
         x = "Condition", y = "Proportion of cells",  fill = "Cell type")
  list(counts = celltype_counts, plot = plot)}

celltype_diff <- list(P1 = run_celltype_diff(P1, "Patient 1"),
                      P2 = run_celltype_diff(P2, "Patient 2"))

CT_P1con <- celltype_diff$P1$plot
CT_P2con <- celltype_diff$P2$plot

### STEP 5: DE genes within cell types ###
run_DE <- function(seurat_obj, celltype)
  {sub <- subset(seurat_obj, subset = manual_celltype == celltype)
  Idents(sub) <- "condition" 
  list(object = sub, AvN = FindMarkers(sub, ident.1 = "Adenoma", ident.2 = "Normal"),
       CvN = FindMarkers(sub, ident.1 = "Carcinoma", ident.2 = "Normal"),
       AvC = FindMarkers(sub, ident.1 = "Adenoma", ident.2 = "Carcinoma"))}

de_results <- list(P1_T = run_DE(P1, "T cells"),
                   P1_epi = run_DE(P1, "epithelial"),
                   P2_T = run_DE(P2, "T cells"),
                   P2_epi = run_DE(P2, "epithelial"))

DE_P1T <- FeaturePlot(de_results$P1_T$object, features = c("CD3D", "CD3E", "TRAC"), 
                      split.by = "condition")
DE_P2T <- FeaturePlot(de_results$P2_T$object, features = c("CD3D", "CD3E", "TRAC"), 
                      split.by = "condition")
DE_P1E <- FeaturePlot(de_results$P1_epi$object, features = c("EPCAM", "KRT8", "KRT18", "KRT19"), 
                      split.by = "condition")
DE_P2E <- FeaturePlot(de_results$P2_epi$object, features = c("EPCAM", "KRT8", "KRT18", "KRT19"), 
                      split.by = "condition")


### STEP 6: clusterProfiler - enrichment analysis ###
run_enrichment <- function(de_table) {de_table$gene <- rownames(de_table)
  sig_genes <- de_table %>%
    dplyr::filter(p_val_adj < 0.01, abs(avg_log2FC) > 1)
  gene_entrez <- bitr(sig_genes$gene, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene = gene_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "hsa", 
                      pvalueCutoff = 0.05)
  list(genes = sig_genes, entrez = gene_entrez, GO = ego, KEGG = ekegg,
    GO_plot = dotplot(ego, showCategory = 15),
    KEGG_plot = dotplot(ekegg, showCategory = 15))}

enrichment_results <- list()
for(sample_name in names(de_results)) {
  enrichment_results[[sample_name]] <- list()
  for(comparison in c("AvN", "CvN", "AvC")) {
    enrichment_results[[sample_name]][[comparison]] <-
      run_enrichment(de_results[[sample_name]][[comparison]])}}

#GO enrichment analysis
plot_GO_comparisons <- function(results, sample_name){
  results[[sample_name]]$AvN$GO_plot +
    results[[sample_name]]$CvN$GO_plot +
    results[[sample_name]]$AvC$GO_plot +
    plot_layout(ncol = 3)}

GO_P1T <- plot_GO_comparisons(enrichment_results, "P1_T")
GO_P2T <- plot_GO_comparisons(enrichment_results, "P2_T")
GO_P1E <- plot_GO_comparisons(enrichment_results, "P1_epi")
GO_P2E <- plot_GO_comparisons(enrichment_results, "P2_epi")

#KEGG enrichment
plot_KEGG_comparisons <- function(results, sample_name){
  results[[sample_name]]$AvN$KEGG_plot +
    results[[sample_name]]$CvN$KEGG_plot +
    results[[sample_name]]$AvC$KEGG_plot +
    plot_layout(ncol = 3)}

KEGG_P1T <-plot_KEGG_comparisons(enrichment_results, "P1_T")
KEGG_P2T <-plot_KEGG_comparisons(enrichment_results, "P2_T")
KEGG_P1E <-plot_KEGG_comparisons(enrichment_results, "P1_epi")
KEGG_P2E <-plot_KEGG_comparisons(enrichment_results, "P2_epi")


### STEP 7: GSEA ###
GSEA_enrichment <- function(de_table) {de_table$gene <- rownames(de_table)
  gene_ranking <- de_table$avg_log2FC
  names(gene_ranking) <- de_table$gene
  gene_ranking <- gene_ranking[!is.na(gene_ranking)]
  gene_ranking <- sort(gene_ranking, decreasing = TRUE)
  gene_df <- bitr(names(gene_ranking), fromType = "SYMBOL", toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)
  ranking_df <- data.frame(SYMBOL = names(gene_ranking), logFC = gene_ranking)
  ranking_df <- ranking_df %>%
    inner_join(gene_df, by = "SYMBOL") %>%
    distinct(ENTREZID, .keep_all = TRUE)
  gene_list <- ranking_df$logFC
  names(gene_list) <- ranking_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  gsea_go <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP", 
                   keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500, 
                   pvalueCutoff = 0.05, verbose = FALSE)
  gsea_kegg <- gseKEGG(geneList = gene_list, organism = "hsa", minGSSize = 10, 
                       maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  gsea_go_sim <- pairwise_termsim(gsea_go)
  gsea_go_sim@result <- gsea_go_sim@result[!grepl("renal", gsea_go_sim@result$Description),]
  list(ranking_df = ranking_df, gene_list = gene_list, gsea_go = gsea_go, 
       gsea_kegg = gsea_kegg,gsea_go_sim = gsea_go_sim,
    GSEA_GO_plot = dotplot(gsea_go, showCategory = 10),
    GSEA_KEGG_plot = dotplot(gsea_kegg, showCategory = 10),
    GSEA_emapplot = emapplot(gsea_go_sim, showCategory = 10))}

GSEA_results <- list()
for(sample_name in names(de_results)) {
  GSEA_results[[sample_name]] <- list()
  for(comparison in c("AvN", "CvN", "AvC")) {
    GSEA_results[[sample_name]][[comparison]] <-
      GSEA_enrichment(de_results[[sample_name]][[comparison]])}}

#GSEA for GO enrichment analysis
plot_GSEA_GO_comparisons <- function(results, sample_name){
  results[[sample_name]]$AvN$GSEA_GO_plot +
    results[[sample_name]]$CvN$GSEA_GO_plot +
    results[[sample_name]]$AvC$GSEA_GO_plot +
    plot_layout(ncol = 3)}

GSEA_GO_P1T <- plot_GSEA_GO_comparisons(GSEA_results, "P1_T")
GSEA_GO_P2T <- plot_GSEA_GO_comparisons(GSEA_results, "P2_T")
GSEA_GO_P1E <- plot_GSEA_GO_comparisons(GSEA_results, "P1_epi")
GSEA_GO_P2E <- plot_GSEA_GO_comparisons(GSEA_results, "P2_epi")

#GSEA for KEGG enrichment:
plot_GSEA_KEGG_comparisons <- function(results, sample_name){
  results[[sample_name]]$AvN$GSEA_KEGG_plot +
    results[[sample_name]]$CvN$GSEA_KEGG_plot +
    results[[sample_name]]$AvC$GSEA_KEGG_plot +
    plot_layout(ncol = 3)}

GSEA_KEGG_P1T <- plot_GSEA_KEGG_comparisons(GSEA_results, "P1_T")
GSEA_KEGG_P2T <- plot_GSEA_KEGG_comparisons(GSEA_results, "P2_T")
GSEA_KEGG_P1E <- plot_GSEA_KEGG_comparisons(GSEA_results, "P1_epi")
GSEA_KEGG_P2E <- plot_GSEA_KEGG_comparisons(GSEA_results, "P2_epi")

# pathway enrichment
plot_GSEA_emapplot_comparisons <- function(results, sample_name){
  results[[sample_name]]$AvN$GSEA_emapplot +
    results[[sample_name]]$CvN$GSEA_emapplot +
    results[[sample_name]]$AvC$GSEA_emapplot +
    plot_layout(ncol = 3)}

GSEA_ema_P1T <- plot_GSEA_emapplot_comparisons(GSEA_results, "P1_T")
GSEA_ema_P2T <- plot_GSEA_emapplot_comparisons(GSEA_results, "P2_T")
GSEA_ema_P1E <-plot_GSEA_emapplot_comparisons(GSEA_results, "P1_epi")
GSEA_ema_P2E <-plot_GSEA_emapplot_comparisons(GSEA_results, "P2_epi")


### STEP 8: Cross-patient pathway reproducibility ###
#T-cells
extract_gsea <- function(gsea_obj, comparison) {gsea_obj@result %>%
    dplyr::select(Description, NES, p.adjust) %>%
    dplyr::mutate(comparison = comparison)}

cross_patient_reproducibility <- function(GSEA_results, celltype) {
  gsea_all <- bind_rows(
    extract_gsea(GSEA_results[[paste0("P1_", celltype)]]$AvN$gsea_go, "P1_AvN"),
    extract_gsea(GSEA_results[[paste0("P1_", celltype)]]$CvN$gsea_go, "P1_CvN"),
    extract_gsea(GSEA_results[[paste0("P1_", celltype)]]$AvC$gsea_go, "P1_AvC"),
    extract_gsea(GSEA_results[[paste0("P2_", celltype)]]$AvN$gsea_go, "P2_AvN"),
    extract_gsea(GSEA_results[[paste0("P2_", celltype)]]$CvN$gsea_go, "P2_CvN"),
    extract_gsea(GSEA_results[[paste0("P2_", celltype)]]$AvC$gsea_go, "P2_AvC"))
  top_terms <- gsea_all %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(max_abs_NES = max(abs(NES), na.rm = TRUE)) %>%
    dplyr::slice_max(max_abs_NES, n = 25) %>%
    dplyr::pull(Description)
  heatmap_df <- gsea_all %>%
    dplyr::filter(Description %in% top_terms) %>%
    dplyr::select(Description, comparison, NES) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = NES)
  heatmap_mat <- as.data.frame(heatmap_df) 
  rownames(heatmap_mat) <- heatmap_mat$Description
  heatmap_mat$Description <- NULL
  heatmap_mat <- as.matrix(heatmap_mat)
  heatmap_mat[is.na(heatmap_mat)] <- 0
  heatmap_mat[is.nan(heatmap_mat)] <- 0
  heatmap_mat[is.infinite(heatmap_mat)] <- 0
  pheatmap(heatmap_mat, cluster_rows = TRUE, cluster_cols = FALSE,
    main = paste("GSEA pathway activity across CRC progression in", celltype),
    na_col = "grey90")}


H_T <- cross_patient_reproducibility(GSEA_results, "T")

H_epi <- cross_patient_reproducibility(GSEA_results, "epi")



###STEP 9: Saving figures ###
#QC
QC_all <- (VP_P1 | FS_P1) / (VP_P2 | FS_P2)
ggsave("figures/QC.png", QC_all + plot_annotation
       (title = "Quality Control Metrics", 
         subtitle = "Patient 1 (top row), Patient 2 (bottom row)"), width = 14, height = 10)


#Normalization
ggsave("figures/Normalization.png", LP_P1 + LP_P2, width = 14, height = 5)

#Clustering + CellType
ggsave("figures/markers_in_clusters.png", TM_P1 + TM_P2, width = 20, height = 8)

ggsave("figures/Cell_Type.png", (UMAP_P1con + UMAP_P2con) / (DP_P1 + DP_P2) 
       + plot_annotation(title = "Clustering + CellType", 
         subtitle = "Patient 1 (left), Patient 2 (right)"), width = 16, height = 12)

ggsave("figures/celltype_diff_con.png", CT_P1con + CT_P2con, width = 20, height = 8)

#GSEA
ggsave("figures/GSEA_GO_T_con.png", GSEA_GO_P1T / GSEA_GO_P2T + 
         plot_annotation(title = "GSEA for T-cell GO enrichment analysis", 
                         subtitle = "Patient 1 (top row), Patient 2 (bottom row)"), 
       width = 15, height = 14)
ggsave("figures/GSEA_GO_Epi_con.png", GSEA_GO_P1E / GSEA_GO_P2E + 
         plot_annotation(title = "GSEA for epithelial GO enrichment analysis", 
                         subtitle = "Patient 1 (top row), Patient 2 (bottom row)"), 
       width = 15, height = 12)

ggsave("figures/GSEA_KEGG_T_con.png", GSEA_KEGG_P1T / GSEA_KEGG_P2T + 
         plot_annotation(title = "GSEA for T-cell KEGG enrichment", 
                         subtitle = "Patient 1 (top row), Patient 2 (bottom row)"), 
       width = 15, height = 14)
ggsave("figures/GSEA_KEGG_Epi_con.png", GSEA_KEGG_P1E / GSEA_KEGG_P2E + 
         plot_annotation(title = "GSEA for epithelial KEGG enrichment", 
                         subtitle = "Patient 1 (top row), Patient 2 (bottom row)"), 
       width = 15, height = 12)

ggsave("figures/GSEA_ema_T_con.png", GSEA_ema_P1T / GSEA_ema_P2T + 
         plot_annotation(title = "GSEA for T-cell pathway enrichment", 
                         subtitle = "Patient 1 (top row), Patient 2 (bottom row)"), 
       width = 15, height = 10)
ggsave("figures/GSEA_ema_Epi_con.png", GSEA_ema_P1E / GSEA_ema_P2E + 
         plot_annotation(title = "GSEA for epithelial pathway enrichment", 
                         subtitle = "Patient 1 (top row), Patient 2 (bottom row)"), 
       width = 15, height = 8)

#Cross patient
ggsave("figures/Cross_patient_T.png", H_T, width = 14, height = 8)
ggsave("figures/Cross_patient_epi.png", H_epi, width = 14, height = 8)



writeLines(capture.output(sessionInfo()),"sessionInfo.txt")

dev.off()
