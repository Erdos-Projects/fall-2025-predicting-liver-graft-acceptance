if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maEndToEnd")

#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)
library(GEOquery)
library(affy)
library(ggplot2)
library(dplyr)



raw_data_dir <- "C:/Users/nanth/Desktop/Erdo_work/GSE145780_RAW"

cel_files <- list.files(raw_data_dir, pattern = "[.]CEL$", full.names = TRUE, ignore.case = TRUE)

rawData <- ReadAffy(filenames = cel_files)

gse <- getGEO("GSE145780", GSEMatrix = FALSE)

# Extract sample (phenotype) data
samples <- GSMList(gse)

pheno <- data.frame(
  sample = names(samples),
  Cluster = sapply(samples, function(x) {
    gsub("cluster \\(1, 2, 3, 4\\): ", "", Meta(x)$characteristics_ch1[2])
  })
)

# Convert Cluster to a factor
pheno$Cluster <- factor(pheno$Cluster)
phenoData(rawData) <- AnnotatedDataFrame(pheno)
sampleNames(protocolData(rawData)) <- sampleNames(rawData)

#Raw PCA
exp_raw <- log2(Biobase::exprs(rawData))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar_raw <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)

dataGG_raw <- data.frame(
  PC1_raw = PCA_raw$x[,1],
  PC2_raw = PCA_raw$x[,2],
  Cluster = Biobase::pData(rawData)$Cluster  # use your Cluster factor
)
  ggplot(dataGG_raw, aes(PC1_raw, PC2_raw)) +
  geom_point(aes(shape = Cluster, colour = Cluster), size = 3) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar_raw[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar_raw[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +     # adjust shapes for number of clusters
  scale_color_manual(values = c("darkorange2", "dodgerblue4", "forestgreen", "purple"))


# Compute RMA, robust multichip average for summarization, without prior normalization
chang_eset <- rma(rawData, target = "core", normalize = FALSE)

# Compute relative log expression by median log2 intensity of every transcript
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(chang_eset)))

RLE_data <- sweep(Biobase::exprs(chang_eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)

RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)



ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 90, size = 6.5, hjust = 1 ,
                                   ))

# RMA calibration: Background-correct, normalize, and summarize
chang_eset_norm <- rma(rawData, target = "core")

# PCA analysis
# 2. Assign phenotype table to your ExpressionSet
Biobase::pData(chang_eset_norm) <- pheno
exp_chang <- Biobase::exprs(chang_eset_norm)
PCA <- prcomp(t(exp_chang), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(
  PC1 = PCA$x[,1],
  PC2 = PCA$x[,2],
  Cluster = Biobase::pData(chang_eset_norm)$Cluster
)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Cluster, colour = Cluster), size = 3) +  # shape and color by Cluster
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +    # one shape per cluster
  scale_color_manual(values = c("darkorange2", "dodgerblue4", "forestgreen", "purple"))  # one color per cluster

# Heatmap clustering, sample-to-sample distances
library(pheatmap)
library(RColorBrewer)
clusters <- Biobase::pData(chang_eset_norm)$Cluster

annotation_for_heatmap <- data.frame(
  Cluster = Biobase::pData(chang_eset_norm)$Cluster
)

row.names(annotation_for_heatmap) <- row.names(Biobase::pData(chang_eset_norm))

dists <- as.matrix(dist(t(exp_chang), method = "manhattan"))

rownames(dists) <- row.names(pData(chang_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

# Annotation colors for Cluster
ann_colors <- list(
  Cluster = c(
    "1" = "darkorange2",
    "2" = "dodgerblue4",
    "3" = "forestgreen",
    "4" = "purple"
  )
)

pheatmap(dists, 
         col = hmcol, 
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), max(dists, na.rm = TRUE)),
         legend_labels = c("small distance", "large distance"),
         main = "Clustering heatmap grouped by Cluster")


# Filtering based on intensity
chang_medians <- rowMedians(Biobase::exprs(chang_eset_norm))

hist_res <- hist(chang_medians, breaks = 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities",
                 xlim = c(min(chang_medians) - 0.5, max(chang_medians) + 0.5))  # expand by 1 unit on each side


     # Visually set cutoff line
      man_threshold <- 3.5
      
      hist_res <- hist(chang_medians, breaks = 100, col = "cornsilk1", freq = FALSE,
                       main = "Histogram of the median intensities",
                       border = "antiquewhite4",
                       xlab = "Median intensities",
                       xlim = c(min(chang_medians) - 0.5, max(chang_medians) + 0.5))  # expand by 1 unit on each side
    
      abline(v = man_threshold, col = "coral4", lwd = 2)
      
# Finish step 10 on filtering by intensity and step 11
no_of_samples <- 
  table(pData(chang_eset_norm)$Cluster)
no_of_samples 

samples_cutoff <- min(no_of_samples)

# Change man_threshold to fit needs
idx_man_threshold <- apply(Biobase::exprs(chang_eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

chang_manfiltered <- subset(chang_eset_norm, idx_man_threshold)

# Annotation
# Load GPL15207 annotation, skipping comment lines
annotation_file <- "C:/Users/nanth/Desktop/Erdo_work/GPL15207-17536.txt"
annotation_data <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#"
)

# Keep only probes with gene symbols
anno_chang <- annotation_data %>%
  filter(Gene.Symbol != "" & !is.na(Gene.Symbol)) %>%
  select(PROBEID = ID, Gene.Symbol, Gene.Title)

anno_chang <- anno_chang %>%
  mutate(Gene.Symbol = sapply(strsplit(Gene.Symbol, " /// "), `[`, 1))

# Merge expression data with annotation
expr_df <- as.data.frame(exprs(chang_manfiltered))
expr_df$PROBEID <- rownames(expr_df)

expr_gene <- merge(expr_df, anno_chang, by.x = "PROBEID", by.y = "PROBEID")

# Collapse multiple probes mapping to the same gene (median)
expr_gene <- expr_gene %>%
  filter(Gene.Symbol != "---" & !is.na(Gene.Symbol))

expr_gene_summary <- expr_gene %>%
  group_by(Gene.Symbol) %>%
  summarize(across(where(is.numeric), median), .groups = "drop")

# Export to CSV with Gene.Symbol as the first column
write.csv(expr_gene_summary, 
          file = "C:/Users/nanth/Desktop/Erdo_work/expr_gene_summary.csv", 
          row.names = FALSE)

# Extract sample metadata
sample_info <- pData(chang_eset_norm)

colnames(sample_info)

sample_cluster_df <- data.frame(
  Sample = rownames(sample_info),
  Cluster = sample_info$Cluster
)

write.csv(sample_cluster_df, 
          file = "C:/Users/nanth/Desktop/Erdo_work/sample_cluster.csv", 
          row.names = FALSE)