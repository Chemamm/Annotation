# Install and load the gplots package if not already installed
library(pheatmap)

# Check if the correct number of command-line arguments is provided
if (length(commandArgs(TRUE)) != 2) {
  cat("Usage: Rscript heatmap_script.R <input_matrix_file> <output_heatmap_file>\n")
  quit(status = 1)
}

# Get input matrix file and output heatmap file from command-line arguments
expression_matrix <- commandArgs(TRUE)[1]
output_heatmap_file <- commandArgs(TRUE)[2]

expression_matrix <- read.delim(expression_matrix, header = TRUE, sep=",", row.names = 1, stringsAsFactors
 = FALSE)



png(paste(output_heatmap_file, "_all_scale.png", sep=""), width=2000, height=2000, res=300)
pheatmap(expression_matrix,
         color = colorRampPalette(c("blue", "grey", "red"))(100),  # Specify the color scheme
         show_rownames = TRUE,                      # Show gene names on the heatmap
         show_colnames = TRUE,                      # Show sample names on the heatmap
         fontsize = 5,                             # Adjust the font size of labels
         scale="row",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
)

dev.off()

png(paste(output_heatmap_file, "_all.png", sep=""), width=2000, height=2000, res=300)
pheatmap(expression_matrix,
         color = colorRampPalette(c("grey", "red"))(100),  # Specify the color scheme
         show_rownames = TRUE,                      # Show gene names on the heatmap
         show_colnames = TRUE,                      # Show sample names on the heatmap
         fontsize = 5,
         cluster_cols=FALSE,
         cluster_rows=FALSE,
)

dev.off()

sg_samples <- colnames(expression_matrix)[grepl("SG_", colnames(expression_matrix))]
mg_samples <- colnames(expression_matrix)[grepl("MG_", colnames(expression_matrix))]

sg_mean <- rowMeans(expression_matrix[, sg_samples])
mg_mean <- rowMeans(expression_matrix[, mg_samples])




png(paste(output_heatmap_file, "_tissue_scale.png", sep=""), width=2000, height=2000, res=300)
pheatmap(cbind(SG = sg_mean, MG = mg_mean),
         color = colorRampPalette(c("blue", "grey", "red"))(100),  # Specify the color scheme
         show_rownames = TRUE,                      # Show gene names on the heatmap
         show_colnames = TRUE,                      # Show sample names on the heatmap
         fontsize = 5,                             # Adjust the font size of labels
         scale="row",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
)

dev.off()

png(paste(output_heatmap_file, "_tissue.png", sep=""), width=2000, height=2000, res=300)
pheatmap(cbind(SG = sg_mean, MG = mg_mean),
         color = colorRampPalette(c("grey", "red"))(100),  # Specify the color scheme
         show_rownames = TRUE,                      # Show gene names on the heatmap
         show_colnames = TRUE,                      # Show sample names on the heatmap
         fontsize = 5,
         cluster_cols=FALSE,
         cluster_rows=FALSE,
)

dev.off()


### logarithmic scale
log_mat <- log10(expression_matrix)

# Replace "-Inf" values in the log-transformed matrix with NaN
log_mat[log_mat == -Inf] <- NaN

    print(log_mat)

png(paste(output_heatmap_file, "_all_log.png", sep=""), width=2000, height=2000, res=300)
pheatmap(log_mat,
         color = colorRampPalette(c("blue","grey","red"))(100),  # Specify the color scheme
         show_rownames = TRUE,                      # Show gene names on the heatmap
         show_colnames = TRUE,                      # Show sample names on the heatmap
         fontsize = 5,
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         na_col ="black",
)

dev.off()

log_sg_mean <- log10(sg_mean)
log_mg_mean <- log10(mg_mean)
log_sg_mean[log_sg_mean == -Inf] <- NaN
log_mg_mean[log_mg_mean == -Inf] <- NaN

png(paste(output_heatmap_file, "_tissue_log.png", sep=""), width=2000, height=2000, res=300)
pheatmap(cbind(SG = log_sg_mean, MG = log_mg_mean),
         color = colorRampPalette(c("blue","grey","red"))(100),  # Specify the color scheme
         show_rownames = TRUE,                      # Show gene names on the heatmap
         show_colnames = TRUE,                      # Show sample names on the heatmap
         fontsize = 5,
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         na_col ="black",
)

dev.off()




