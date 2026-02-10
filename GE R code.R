library(reshape2)
library(limma)
library(cluster) 
library(gplots)
library(RColorBrewer)
library(matrixStats)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
#library(ggplot2)



# Read the data files
annotations <- readRDS("STA5MB_2025_BC_annotations.RDS")
clinical <- readRDS("STA5MB_2025_BC_clinical_data.RDS")
expression <- readRDS("STA5MB_2025_BC_expression_data.RDS")

# Basic data exploration
cat("annotations class:", class(annotations), "\n")
cat("clinical class:", class(clinical), "\n")
cat("expr_raw class:", class(expression), "\n\n")

cat("Dataset Dimensions:\n")
cat("Annotations:", dim(annotations), "\n")
cat("Clinical data:", dim(clinical), "\n")
cat("Expression data:", dim(expression), "\n")

# Check each dataset structure
cat("\nData Structure:\n")
str(annotations)
str(clinical)
#str(expression)

summary(annotations)
#is.na(annotations)
summary(clinical)
#summary(expression)
head(expression)

cat("# of missing value in Expression:", sum(rowSums(is.na(expression))), "\n")
cat("# of missing value in Clinical:", sum(colSums(is.na(clinical))), "\n")
# 30 Null value, Prepare to clean for Survival Analysis 
cat("# of missing value in Annotation:", sum(colSums(is.na(annotations))), "\n")

# Check clinical data & expression match with "sampleID"
setequal(clinical$sampleID, colnames(expression))

# Re-order clinical data to match with the order in expression 
clinical <- clinical[match(colnames(expression), clinical$sampleID),]
# common_ID <- intersect(clinical$sampleID, colnames(expression))

boxplot(expression, outline = FALSE, main = "Expression distribution per sample",
        ylab = "Expression", xlab = "", las = 2, col = "lightblue")

#hist(as.numeric(expression), breaks = 100, main = "Overall Expression Value Distribution",
#     xlab = "Expression (log2)", col = "skyblue")


expr_norm <- normalizeBetweenArrays(expression, method = "scale") 
boxplot(expr_norm, outline = FALSE, main = "Expression distribution after Normalization",
        ylab = "Expression", xlab = "", las = 2, col = "blue")

#------------------------------------------
# plot avg gene expression distribution 
#------------------------------------------

# Calculate gene expression statistics
gene_means <- rowMeans(expr_norm)
gene_medians <- apply(expr_norm, 1, median)

# Plot the distribution
par(mfrow = c(1,3 ))

# Histogram of mean expression
hist(gene_means, breaks = 50, main = "Distribution of Gene Means", 
     xlab = "Mean Expression", col = "lightblue")
abline(v = quantile(gene_means, c(0.10)), col = "red", lty = 2)
# 10% is about 4, which is reasonable cut off line for low expr

cat("Initial dimensions:", dim(expr_norm), "\n")

#-----------------------------------
# Remove control probes 
#-----------------------------------
control_probes <- grep("^AFFX", rownames(expr_norm))
if(length(control_probes) > 0) {
  expr_filtered <- expr_norm[-control_probes, ]
  cat("Removed", length(control_probes), "control probes\n")
} else {
  expr_filtered <- expr_norm
  cat("No AFFX control probes found\n")
}
gene_filter_means <- rowMeans(expr_filtered)

# Histogram of mean expression after remove control probes 
hist(gene_filter_means, breaks = 50, main = "After control probe", 
     xlab = "Mean Expression", col = "blue")
abline(v = quantile(gene_filter_means, c(0.10)), col = "red", lty = 2)

cat("Dimensions after control probes removed:", dim(expr_filtered), "\n")

#--------------------------------
# Remove lowe expr gene
#--------------------------------
#  Using variance 
#gene_variance <- apply(expr_filtered, 1, var)
gene_variance <- rowVars(expr_filtered)
variance_threshold <- quantile(gene_variance, probs = 0.10)  # Bottom 10%

# Using Median Absolute Deviation 
gene_mad <- rowMads(expr_filtered) #apply(expr_filtered, 1, mad)  # MAD = median absolute deviation
mad_threshold <- quantile(gene_mad, probs = 0.10)  # Bottom 10%
#gene_means <- rowMeans(expr_filtered)
#mean_threshold <- quantile(gene_means, 0.20)  # Remove bottom 20%
#gene_means > mean_threshold
normal_genes <- gene_variance > variance_threshold | gene_mad > mad_threshold
expr_final <- expr_filtered[normal_genes, ]

gene_final_means <- rowMeans(expr_final)
# Histogram of mean expression after remove low expr
hist(gene_final_means, breaks = 50, main = "After low expre gene", 
     xlab = "Mean Expression", col = "darkblue")
abline(v = quantile(gene_final_means, c(0.10)), col = "red", lty = 2)

cat("Dimensions after control probes removed:", dim(expr_final), "\n")

# ------------------------------------------
# Clustering 
# ------------------------------------------

# Scale expression data 
scaled.E <- t(scale(t(expr_final), center = TRUE)) 
# calculate distance
dist <- dist(t(scaled.E))
# **** complete ****
# Perform hierarchical clustering with complete linkage.
hc <- hclust(dist, method = "complete")
K<-2:6
sh<-NULL
for(i in K) {
  sh<-c(sh,median(silhouette(cutree(hc,k=i),dist=dist)[,3],na.rm=T))
}
par(mfrow = c(1,1 ))
#Plotsilhouette
plot(K,sh,type="l",main="Mediansilhouette-complete",xlab="Numberofclusters")

# **** average ****
# Perform hierarchical clustering with average linkage.
hc_avg <- hclust(dist, method = "average")
K<-2:6
sh<-NULL
for(i in K) {
  sh<-c(sh,median(silhouette(cutree(hc_avg,k=i),dist=dist)[,3],na.rm=T))
}
#Plotsilhouette
plot(K,sh,type="l",main="Mediansilhouette-average",xlab="Numberofclusters")

# **** ward.D ****
# Perform hierarchical clustering with ward.D linkage.
hc_w <- hclust(dist, method = "ward.D")
K<-2:6
sh<-NULL
for(i in K) {
  sh<-c(sh,median(silhouette(cutree(hc_w,k=i),dist=dist)[,3],na.rm=T))
}
#Plotsilhouette
plot(K,sh,type="l",main="Mediansilhouette-ward.D",xlab="Numberofclusters")

# **** ward.D2 ****
# Perform hierarchical clustering with ward.D2 linkage.
hc_w2 <- hclust(dist, method = "ward.D2")
K<-2:6
sh<-NULL
for(i in K) {
  sh<-c(sh,median(silhouette(cutree(hc_w2,k=i),dist=dist)[,3],na.rm=T))
}
#Plotsilhouette
plot(K,sh,type="l",main="Mediansilhouette-ward.D2",xlab="Numberofclusters")

#------------------------
# ploting dendrogram 
#------------------------
plot(hc, main = "Hierarchical Clustering of Genes - Complete",
     xlab = "Genes", ylab = "Distance", cex = 0.3, hang = -1)
# Add cutoff line
abline(h = 300, col = "red", lwd = 2, )

plot(hc_avg, main = "Hierarchical Clustering of Genes - Avgerage",
     xlab = "Genes", ylab = "Distance", cex = 0.3, hang = -1)
# Add cutoff line
abline(h = 245, col = "red", lwd = 2, )

plot(hc_w, main = "Hierarchical Clustering of Genes - Ward.D",
     xlab = "Genes", ylab = "Distance", cex = 0.3, hang = -1)
# Add cutoff line
abline(h = 1000, col = "red", lwd = 2, )

plot(hc_w2, main = "Hierarchical Clustering of Genes - Ward.D2",
     xlab = "Genes", ylab = "Distance", cex = 0.3, hang = -1)
# Add cutoff line
abline(h = 550, col = "red", lwd = 2, )

#select optimal clusters for ward.D
#cl<-cutree(hc_w,k=K[1])
# ward.D has the highest distance 
#select optimal clusters for Average
cl<-cutree(hc_w,k=K[1])
# ward.D has the highest distance 
table(cl)
# ------------------------------------------
# create Heatmap with high varaince gene
# ------------------------------------------
# Select high variance genes 
rv <- rowVars(expr_final)
idx <- order(-rv)[1:1000] 
# Specify colour palette 
cols <- brewer.pal(length(unique(cl)), "Set1") 
# Specify colour palette 
#cols <- colors()[seq(8, length(colors()), len = length(unique(cl)))] 

# Produce heatmap 
heatmap.2(scaled.E[idx, ], labCol = cl, trace = "none", ColSideColors = cols[cl],  
          col = redgreen(100), Colv = as.dendrogram(hc_w)) 


# ------------------------------------------
# PCA Visualization for 2 Clusters
# ------------------------------------------ 

idx <- order(-rv)[1:1000] 
# Specify colour palette 
cols <- brewer.pal(length(unique(cl)), "Set1")   

# Extract the data for PCA (use the SCALED data for visualization)
par(bg = "white") 
pc <- princomp(scaled.E[idx, ]) 
summary(pc)

# Or extract specific values:
variance_explained <- pc$sdev^2 / sum(pc$sdev^2) * 100
cat("PC1 explains:", round(variance_explained[1], 2), "%\n")
cat("PC2 explains:", round(variance_explained[2], 2), "%\n")
plot(pc$load[, 1:2], col = cl) 
title("PCs 1 and 2 of cancer data \n coloured by clusters") 

# Create MDS plot, colour coded by cluster 
plotMDS(scaled.E[idx, ], col = as.numeric(cl)) 

# ------------------------------------------
# Gene Expression Analysis 
# ------------------------------------------

#Create design matrix based on clusters
design <- model.matrix(~ factor(cl))

# Fit linear model
fit <- lmFit(expr_final, design)
fit <- eBayes(fit)
summary(fit)
#Testing of cluster difference
de_genes_raw <- topTable(fit, coef = 2, number = Inf, adjust.method = "none")
head(de_genes_raw)
de_genes <- topTable(fit, coef = 2, number = Inf, adjust.method = "fdr")
head(de_genes)

bio_sig_genes_raw <- de_genes_raw[abs(de_genes_raw$logFC) > 1 & de_genes_raw$P.Val < 0.05, ]
nrow(bio_sig_genes_raw)
bio_sig_genes <- de_genes[abs(de_genes$logFC) > 1 & de_genes$adj.P.Val < 0.05, ]
nrow(bio_sig_genes)

# ------------------------------------------
# p-value Histogram (before & after FRD)
# ------------------------------------------
par(mfrow = c(1, 2))

# Raw P-values histogram
hist(de_genes_raw$P.Value, breaks = 50, 
     main = "Raw P-value Distribution",
     xlab = "Raw P-value", col = "lightblue",
     xlim = c(0, 1))
abline(v = 0.05, col = "red", lwd = 2, lty = 2)


# Adjusted P-values (FDR) histogram
hist(de_genes$adj.P.Val, breaks = 50,
     main = "FDR-adjusted P-value Distribution",
     xlab = "FDR-adjusted P-value", col = "lightgreen",
     xlim = c(0, 1))
abline(v = 0.05, col = "red", lwd = 2, lty = 2)


# ------------------------------------------
# Volcano Plot
# ------------------------------------------
par(mfrow = c(1,1))
# Smaller fold change
cat("Raw p < 0.05", 
    sum(de_genes_raw$P.Value < 0.05 ), "\n")
cat("FDR < 0.05", 
    sum(de_genes$adj.P.Val < 0.05 ), "\n")
cat("=== DE ANALYSIS ASSESSMENT ===\n")
cat("Total genes analyzed:", nrow(de_genes), "\n")
cat("Significant DE genes (FDR < 0.05):", sum(de_genes$adj.P.Val < 0.05), "\n")
cat("Significant DE genes (FDR < 0.01):", sum(de_genes$adj.P.Val < 0.01), "\n")
cat("Biologically significantly Up-regulated in cluster:", sum(de_genes$adj.P.Val < 0.05 & de_genes$logFC > 1), "\n")
cat("Biologically significantly Down-regulated in cluster:", sum(de_genes$adj.P.Val < 0.05 & de_genes$logFC < -1), "\n")

# Volcano plot to visualize DE results
volcano_data <- de_genes

plot(volcano_data$logFC, -log10(volcano_data$P.Value),
     xlab = "Log2 Fold Change", ylab = "-log10(Raw P-value)",
     main = "Volcano Plot: Differential Expression",
     pch = 16, cex = 0.6,
     col = ifelse(volcano_data$adj.P.Val < 0.05 & abs(volcano_data$logFC) > 1, 
                  "red", "gray"),
     xlim = c(-max(abs(volcano_data$logFC)), max(abs(volcano_data$logFC))))



# ------------------------------------------
# GO Enrichment
# ------------------------------------------
# Create probe_id column
de_genes$probe_id <- rownames(de_genes)

# Merge DE table w annotation
de_annotated <- merge(de_genes, annotations, by.x = "probe_id", 
                      by.y = "affy_hg_u133_plus_2", all.x = TRUE)

#cat("Total NA vlaues: ", colSums(is.na(de_annotated)))
# Get significant DE genes from cancer analysis
sig_genes <- de_annotated[de_annotated$adj.P.Val < 0.05 & abs(de_annotated$logFC) > 1, ]
write.csv(sig_genes, "sig_genes_results.csv", row.names = TRUE)
cat("Total NA vlaues: ", colSums(is.na(sig_genes)))
cat("Total # of Significant Genes: ", nrow(sig_genes))
# Clean symbol data
sig_genes_symbol <- na.omit(unique(sig_genes$hgnc_symbol))
# Covert hgnc_symbol to ENTREZID
sig_entrez <- bitr(sig_genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis for Biological Process (BP)
go_bp <- enrichGO(gene         = sig_entrez$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "BP",          # Can also be "MF" or "CC"
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)          

head(go_bp)
write.csv(go_bp, "GOenrichment_results.csv", row.names = TRUE)


# ------------------------------------------
# Visualization of GO enrichment
# ------------------------------------------

# Dotplot for Biological Process
dotplot(go_bp, showCategory=15, title="GO Enrichment - Biological Process")

# GO Network Plot (BP example)
cnetplot(go_bp, showCategory=10, foldChange=sig_genes$logFC)

par(mfrow = c(1, 1))

# ------------------------------------------
# Survival Analysis 
# ------------------------------------------
setequal(clinical$sampleID, names(cl))
clinical_clean <- clinical[clinical$sampleID %in% names(cl), ]
clinical_clean$cluster <- cl[match(clinical_clean$sampleID, names(cl))]
cat("Missing Value for Clinical Data:", colSums(is.na(clinical_clean)), "\n")

# Clincial data remove imcompleted record 
clinical_clean <- clinical_clean[!is.na(clinical_clean$Surv_time) & !is.na(clinical_clean$event),]
nrow(clinical_clean)

# ------------------------------------------
# Identify Zero Event Columns and Merge
# ------------------------------------------
unique_values_list <- list()
for (col_name in names(clinical_clean)) {
  unique_values_list[[col_name]] <- unique(clinical_clean[[col_name]])
}
print(unique_values_list)
# ERstatus
er_table <- table(clinical_clean$ERstatus, clinical_clean$event)
cat("\nERstatus event distribution:\n")
print(er_table)

# LNstatus  
ln_table <- table(clinical_clean$LNstatus, clinical_clean$event)
cat("\nLNstatus event distribution:\n")
print(ln_table)
#LN-   127   22
#LN?     9    0
#LN+    45   33
cat("ER? represents", round(4/nrow(clinical_clean)*100, 1), "% of total samples\n")
cat("LN? represents", round(9/nrow(clinical_clean)*100, 1), "% of total samples\n")
#FALSE TRUE
#ER-    25    6
#ER?     4    0
#ER+   152   49
# ------------------------------------------
# Collapse Problematic Categories (LN?-> LN-, ER? -> ER-)

clinical_clean$LNstatus[clinical_clean$LNstatus == "LN?"] <- "LN-"
clinical_clean$LNstatus <- factor(clinical_clean$LNstatus)
clinical_clean$ERstatus[clinical_clean$ERstatus == "ER?"] <- "ER-"
clinical_clean$ERstatus <- factor(clinical_clean$ERstatus)


# ------------------------------------------
# Log Rank Test & K-M Survival Plot
# ------------------------------------------
surv_obj <- Surv(time = clinical_clean$Surv_time, event = clinical_clean$event)
km_comp <- survdiff(surv_obj ~ clinical_clean$cluster) # stratify by gender
km_comp

fit_km <- survfit(surv_obj ~ clinical_clean$cluster)
autoplot(fit_km) + 
  labs(x = "\n Survival Time (Years) since diagnosis ", y = "Survival Probabilities \n", 
       title = "K-M Survival Curves for Breast Cancer Patients \n")

fit_km_ln <- survfit(surv_obj ~ clinical_clean$LNstatus)
summary(fit_km_ln)

autoplot(fit_km_ln) + 
  labs(x = "\n Survival Time (Years) since diagnosis ", y = "Survival Probabilities \n", 
       title = "K-M Survival Curves for Breast Cancer Patients \n")



# ------------------------------------------
# Univariate Cox
# ------------------------------------------
#surv_obj_m2 <- Surv(time = clinical_clean$Surv_time, event = clinical_clean$event)
surv_obj <- Surv(time = clinical_clean$Surv_time, event = clinical_clean$event)
# List of variables to test
variables <- c("cluster", "histgrade", "ERstatus", "PRstatus", "age", "tumor_size_mm", "LNstatus")

# Create a function for univariate analysis
perform_univariate_analysis <- function(data, variables) {
  results <- list()
  
  for(var in variables) {
    # Create formula
    if(is.numeric(data[[var]])) {
      # For continuous variables
      formula <- as.formula(paste("surv_obj ~", var))
    } else {
      # For categorical variables
      formula <- as.formula(paste("surv_obj ~ factor(", var, ")"))
    }
    
    # Fit Cox model
    cox_uni <- coxph(formula, data = data)
    cat("\nCox Result for", var, ":\n")
    # Store results
    print(summary(cox_uni))
  }
}

# Perform univariate analysis (merge cat)
uni_results_2 <- perform_univariate_analysis(clinical_clean, variables)

# ------------------------------------------
# Multivariate Model
# ------------------------------------------

# Adjusted for prognostic factors: histological grade, ER, PR, age, tumor size, lymph nodes
cox_model_c1<- coxph(surv_obj ~ cluster + histgrade + ERstatus + PRstatus + age + tumor_size_mm + LNstatus,
                    data = clinical_clean)

summary(cox_model_c1)

cox_model_c2<- coxph(surv_obj ~ cluster + as.factor(tumor_size_mm) + LNstatus,
                    data = clinical_clean)

summary(cox_model_c2)


