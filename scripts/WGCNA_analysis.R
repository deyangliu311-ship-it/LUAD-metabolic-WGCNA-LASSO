rm(list=ls())
set.seed(123)

############################
# 1. Load libraries
############################
library(TCGAbiolinks)
library(edgeR)
library(limma)
library(tidyverse)
library(GSVA)
library(ConsensusClusterPlus)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(glmnet)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

############################
# 2. Download TCGA LUAD data
############################
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

GDCdownload(query)
data <- GDCprepare(query)

expr_raw <- assay(data)
clinical <- colData(data)

############################
# 3. TPM normalization
############################
dge <- DGEList(counts = expr_raw)
dge <- calcNormFactors(dge)
tpm <- cpm(dge, normalized.lib.sizes = TRUE)
expr <- log2(tpm + 1)

############################
# 4. Metabolic gene selection
############################
metabolic_genes <- read.table("metabolic_genes.txt")$V1
expr_meta <- expr[intersect(rownames(expr), metabolic_genes), ]

############################
# 5. Consensus clustering
############################
cc <- ConsensusClusterPlus(as.matrix(expr_meta),
                           maxK=10, reps=1000,
                           pItem=0.8, pFeature=1,
                           clusterAlg="hc", distance="pearson",
                           seed=1234)

clusters <- cc[[3]]$consensusClass
names(clusters) <- colnames(expr_meta)

############################
# 6. PCA
############################
pca <- prcomp(t(expr_meta), scale.=TRUE)
pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Cluster=as.factor(clusters))

write.csv(pca_df, "PCA_results.csv")

############################
# 7. GSVA pathway scoring
############################
gmt <- getGmt("metabolism_pathways.gmt")
gsva_scores <- gsva(expr, gmt, method="gsva")

write.csv(gsva_scores, "GSVA_scores.csv")

############################
# 8. WGCNA
############################
datExpr <- t(expr)
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector=powers)
softPower <- 5

net <- blockwiseModules(datExpr, power=softPower,
                        TOMType="unsigned",
                        minModuleSize=30,
                        mergeCutHeight=0.25,
                        numericLabels=TRUE,
                        verbose=3)

moduleColors <- labels2colors(net$colors)

############################
# 9. Hub genes (red module)
############################
hub_genes <- names(moduleColors[moduleColors=="red"])

############################
# 10. Enrichment analysis
############################
ego <- enrichGO(hub_genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL")
ekegg <- enrichKEGG(hub_genes, organism="hsa")

write.csv(as.data.frame(ego), "GO_enrichment.csv")
write.csv(as.data.frame(ekegg), "KEGG_enrichment.csv")

############################
# 11. Survival data preparation
############################
clinical$OS.time <- as.numeric(clinical$days_to_death)
clinical$OS.status <- ifelse(clinical$vital_status=="Dead",1,0)

common_samples <- intersect(colnames(expr), clinical$barcode)
expr_surv <- expr[hub_genes, common_samples]

surv_df <- data.frame(
  time = clinical$OS.time[match(common_samples, clinical$barcode)],
  status = clinical$OS.status[match(common_samples, clinical$barcode)],
  t(expr_surv)
)

############################
# 12. LASSO Cox
############################
x <- as.matrix(surv_df[,-c(1,2)])
y <- Surv(surv_df$time, surv_df$status)

lasso <- cv.glmnet(x, y, family="cox")
coef_lasso <- coef(lasso, s="lambda.min")
sel_genes <- rownames(coef_lasso)[coef_lasso!=0]

write.csv(sel_genes, "selected_genes.csv")

############################
# 13. Cox model
############################
cox_df <- surv_df[,c("time","status",sel_genes)]
cox_model <- coxph(Surv(time,status)~., data=cox_df)

risk_score <- predict(cox_model, type="risk")
group <- ifelse(risk_score>median(risk_score),"High","Low")

############################
# 14. Kaplan-Meier
############################
fit <- survfit(Surv(time,status)~group)

pdf("KM_curve.pdf")
ggsurvplot(fit, data=data.frame(group))
dev.off()

############################
# 15. ROC curves
############################
roc1 <- timeROC(surv_df$time, surv_df$status, risk_score, times=365)
roc3 <- timeROC(surv_df$time, surv_df$status, risk_score, times=1095)
roc5 <- timeROC(surv_df$time, surv_df$status, risk_score, times=1825)

pdf("ROC_curves.pdf")
plot(roc1$FP, roc1$TP, col="red", type="l")
lines(roc3$FP, roc3$TP, col="blue")
lines(roc5$FP, roc5$TP, col="green")
legend("bottomright", legend=c("1-year","3-year","5-year"),
       col=c("red","blue","green"), lty=1)
dev.off()

############################
# 16. Heatmap
############################
pdf("heatmap_genes.pdf")
pheatmap(expr[sel_genes, common_samples],
         annotation_col=data.frame(Risk=group))
dev.off()

############################
# END
############################
