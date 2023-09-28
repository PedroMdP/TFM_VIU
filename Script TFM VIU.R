library(Biobase)
library(limma)
library(ggVennDiagram)
library(ggplot2)
library(gplots)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(cowplot)


#Establecimiento del directorio de trabajo
setwd("/Users/pedro/Documents/VIU/0 - Trabajo Fin de Máster/Experimentación/Lecturas prueba")

#Lectura metadatos, creación de grupos y tabla "diseño"
metadata <- read.delim("metadata.csv", sep = ";", row.name = 1)
group <- paste(metadata$diagnosis, sep = ";")
group <- factor(group)
table(group)
design <- model.matrix(~0+group)
design
colnames(design) <- levels(group)
rownames(design) <- rownames(metadata)
diagnostic <- paste(metadata$diagnosis)
diagnostic

#Establecimiento del contraste
contrast <- makeContrasts(SS-NSS, levels = group)
contrast

#Anotación ID y nombre gen
featuredata <- read.table("GPL13497-9755.txt", sep = "\t", header = TRUE, fill = TRUE, quote = "\"")
featuredata <- featuredata[, c(1,7)]
featuredata$ID2 <- featuredata$ID
featuredata <- column_to_rownames(featuredata, var = "ID2")

#Sustracción de las sondas
featuredata <- featuredata[grep("A_2|A_3", featuredata$ID) , ]

#Creación de lista de lectura de archivos .txt de microarrays
targets <- readTargets("Targets.txt", row.names = "Filename")

#Lectura archivos de arrays
x <- read.maimages(targets$Filename, source = "agilent.median", green.only = TRUE)

#Control de calidad de las lecturas antes de la eliminación de fondo y de la normalización 
plot_no_norm <- plotDensities(x, col = "black", legend = FALSE)
plotMD1_no_norm <- plotMD(x, column = 1) + abline(h = 0, col = "red", lty = 2, lwd = 2)
plotMD2_no_norm <- plotMD(x, column = 2) + abline(h = 0, col = "red", lty = 2, lwd = 2)

#Corrección de fondo y normalización
x_fondo <- backgroundCorrect(x, method = "normexp")



#Método de normalización "none"

x_norm_none <- normalizeBetweenArrays(x_fondo, method = "none")

#Control distribución después de la eliminación de fondo y de la normalización 
plot_norm_none <- plotDensities(x_norm_none, col = "black", legend = FALSE)
plotMD1_norm_none <- plotMD(x_norm_none, column = 1) + abline(h = 0, col="red", lty = 2, lwd = 2)
plotMD2_norm_none <- plotMD(x_norm_none, column = 2) + abline(h = 0, col="red", lty = 2, lwd = 2)

#Transformación en base de datos
data_none <- as.data.frame(x_norm_none)

#Sustracción de las sondas
data_none <- data_none[grep("A_2|A_3", data_none$ProbeName) , ]
data_none <- data_none[ ,-c(1:3,5)]

#Agregación de la base de datos por método de la media en función de ProbeName
data_agreg_none <- aggregate(data_none, by = list(data_none$ProbeName), FUN = mean, na.rm = TRUE)
data_agreg_none <- data_agreg_none[ ,-c(2)]
colnames(data_agreg_none)[1] = "ID"

#Anotación de la base de datos
data_anot_none <- data_agreg_none %>% inner_join(., featuredata, by = "ID")

#Eliminación de las sondas que no tienen correlación con ningún gen
data_anot_none_noNA <- data_anot_none[!data_anot_none$GENE_SYMBOL == "",]

#Segunda agregación de la base de datos por método de la media en función de GENE_SYMBOL
data_final_none <- aggregate(data_anot_none_noNA[, -c(1,131)], by = list(data_anot_none_noNA$GENE_SYMBOL), FUN = mean, na.rm = TRUE)
colnames(data_final_none)[1] = "Gene"
rownames(data_final_none) <- NULL
data_final_none <- column_to_rownames(data_final_none, var = "Gene")

#Eliminación de las filas que contienen algún valor negativo
data_final_none[data_final_none < 0] <- NA
data_final_none <- na.omit(data_final_none)

#Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_none <- prcomp(t(data_final_none), scale = TRUE)
plot(pca_none$x[, 1], pca_none$x[, 2])
pca.var_none <- pca_none$sdev^2
pca.var.per_none <- round(pca.var_none/sum(pca.var_none)*100, 1)
barplot(pca.var.per_none, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_none <- data.frame(Sample = rownames(pca_none$x), X = pca_none$x[, 1], Y = pca_none$x[, 2])
ggplot(data = pca.data_none, aes(x = pca_none$x[, 1], y = pca_none$x[, 2], label = Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep = "")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_none <- pca_none$rotation[, 1]
gene_score_none <- abs(loading_score_none)
gene_score_ranked_none <- sort(gene_score_none, decreasing=TRUE)
top_10_genes_none <- names(gene_score_ranked_none[1:10])
top_10_genes_none
pca_none$rotation[top_10_genes_none, 1]
pca.data_none2 <- data.frame(Sample = rownames(pca_none$x), X = pca_none$x[,1], Y = pca_none$x[, 2], Condition = diagnostic)
PCA_none <- ggplot(data = pca.data_none2, aes(x = pca_none$x[, 1], y = pca_none$x[, 2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep = "")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep = "")) + theme_bw() + ggtitle("NONE")
PCA_none

#limma
fitlimma_none <- lmFit(data_final_none, design)
fitlimma_none <- contrasts.fit(fitlimma_none, contrast)
fitlimma_none <- eBayes(fitlimma_none)
toptablelimma_none <- topTable(fitlimma_none, sort.by = "logFC", number = Inf, adjust.method = "none", coef = 1)
toptablelimma_none
results_none <- decideTests(fitlimma_none)

#Volcano plot
# Si el gen cumple logFC > 0.5 y P.Value < 0.05 se denomina como "UP" 
toptablelimma_none$diffexpressed[toptablelimma_none$logFC > 0.5 & toptablelimma_none$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_none$diffexpressed[toptablelimma_none$logFC < -0.5 & toptablelimma_none$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_none$diffexpressed[is.na(toptablelimma_none$diffexpressed)] <- "NO_none"
#Función gráfica Volcano
vp_none <- ggplot(data=toptablelimma_none, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col = "black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("NONE") + theme(plot.title = element_text(hjust = 0.5))
vp_none

#Mapa de calor
#Selección de los 20 primeros genes que mayor diferencia de expresión presentan
names_none <- rownames(toptablelimma_none[1:20, ])
df_none_short <- data_final_none[c(names_none),c(1:114)]
#Anotación de las muestras
annotation_df <- metadata[1:114, 4, drop = FALSE]
#Función mapa de calor
heatmap_none <- pheatmap(df_none_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "None", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row")



#Método de normalización "quantile"

x_norm_quantile <- normalizeBetweenArrays(x_fondo, method = "quantile")

#Control distribución después de normalizar y eliminar fondo 
plot_norm_quantile <- plotDensities(x_norm_quantile, col = "black", legend = FALSE)
plotMD1_norm_quantile <- plotMD(x_norm_quantile, column = 1) + abline(h = 0, col = "red", lty = 2, lwd = 2)
plotMD2_norm_quantile <- plotMD(x_norm_quantile, column = 2) + abline(h = 0, col = "red", lty = 2, lwd = 2)

#Transformación en base de datos
data_quantile <- as.data.frame(x_norm_quantile)

#Sustracción de las sondas
data_quantile <- data_quantile[grep("A_2|A_3", data_quantile$ProbeName) ,]
data_quantile <- data_quantile[ ,-c(1:3,5)]

#Agregación de la base de datos por método de la media en función de ProbeName
data_agreg_quantile <- aggregate(data_quantile, by = list(data_quantile$ProbeName), FUN = mean, na.rm = TRUE)
data_agreg_quantile <- data_agreg_quantile[ ,-c(2)]
colnames(data_agreg_quantile)[1] = "ID"

#Anotación de la base de datos
data_anot_quantile <- data_agreg_quantile %>% inner_join(., featuredata, by = "ID")

#Eliminación de las sondas que no tienen correlación con ningún gen
data_anot_quantile_noNA <- data_anot_quantile[!data_anot_quantile$GENE_SYMBOL=="",]

#Segunda agregación de la base de datos por método de la media en función de GENE_SYMBOL
data_final_quantile <- aggregate(data_anot_quantile_noNA[, -c(1,131)], by = list(data_anot_quantile_noNA$GENE_SYMBOL), FUN = mean, na.rm = TRUE)
colnames(data_final_quantile)[1] = "Gene"
rownames(data_final_quantile) <- NULL
data_final_quantile <- column_to_rownames(data_final_quantile, var = "Gene")

#Eliminación de las filas que contienen algún valor negativo
data_final_quantile[data_final_quantile < 0] <- NA
data_final_quantile <- na.omit(data_final_quantile)

#Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_quantile <- prcomp(t(data_final_quantile), scale=TRUE)
plot(pca_quantile$x[, 1], pca_quantile$x[, 2])
pca.var_quantile <- pca_quantile$sdev^2
pca.var.per_quantile <- round(pca.var_quantile/sum(pca.var_quantile)*100, 1)
barplot(pca.var.per_quantile, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_quantile <- data.frame(Sample = rownames(pca_quantile$x), X = pca_quantile$x[, 1], Y = pca_quantile$x[, 2])
ggplot(data = pca.data_quantile, aes(x = pca_quantile$x[, 1], y = pca_quantile$x[, 2], label = Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_quantile[1], "%", sep = "")) + ylab(paste("PC2 - ", pca.var.per_quantile[2], "%", sep = "")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_quantile <- pca_quantile$rotation[, 1]
gene_score_quantile <- abs(loading_score_quantile)
gene_score_ranked_quantile <- sort(gene_score_quantile, decreasing = TRUE)
top_10_genes_quantile <- names(gene_score_ranked_quantile[1:10])
top_10_genes_quantile
pca_quantile$rotation[top_10_genes_quantile, 1]
pca.data2_quantile <- data.frame(Sample = rownames(pca_quantile$x), X = pca_quantile$x[, 1], Y = pca_quantile$x[, 2], Condition = diagnostic)
PCA_quantile <- ggplot(data = pca.data2_quantile, aes(x = pca_quantile$x[,1], y = pca_quantile$x[,2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_quantile[1], "%", sep = "")) + ylab(paste("PC2 - ", pca.var.per_quantile[2], "%", sep = "")) + theme_bw() + ggtitle("QUANTILE")
PCA_quantile

#limma
fitlimma_quantile <- lmFit(data_final_quantile, design)
fitlimma_quantile <- contrasts.fit(fitlimma_quantile, contrast)
fitlimma_quantile <- eBayes(fitlimma_quantile)
toptablelimma_quantile <- topTable(fitlimma_quantile, sort.by = "logFC", number = Inf, adjust.method = "none", coef = 1)
toptablelimma_quantile
results_quantile <- decideTests(fitlimma_quantile)

#Volcano plot
# Si el gen cumple logFC > 0.5 y  P.Value < 0.05 se denomina como "UP" 
toptablelimma_quantile$diffexpressed[toptablelimma_quantile$logFC > 0.5 & toptablelimma_quantile$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_quantile$diffexpressed[toptablelimma_quantile$logFC < -0.5 & toptablelimma_quantile$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_quantile$diffexpressed[is.na(toptablelimma_quantile$diffexpressed)] <- "NO_none"
#Función gráfica Volcano
vp_quantile <- ggplot(data=toptablelimma_quantile, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col = "black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("QUANTILE") + theme(plot.title = element_text(hjust = 0.5))
vp_quantile

#mapa de calor
#Selección de las 20 primeros genes que mayor diferencia de expresión presentan
names_quantile <- rownames(toptablelimma_quantile[1:20, ])
df_quantile_short <- data_final_quantile[c(names_quantile), c(1:114)]
#Función del mapa de calor
heatmap_quantile <- pheatmap(df_quantile_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Quantile", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row")



#MÉTODO NORMALIZACIÓN "SCALE"

x_norm_scale <- normalizeBetweenArrays(x_fondo, method = "scale")

#Control distribución después de normalizar y eliminar fondo 
plot_norm_scale <- plotDensities(x_norm_scale, col = "black", legend = FALSE)
plotMD1_norm_scale <- plotMD(x_norm_scale, column = 1) + abline(h = 0, col = "red", lty = 2, lwd = 2)
plotMD2_norm_scale <- plotMD(x_norm_scale, column = 2) + abline(h = 0, col = "red", lty = 2, lwd = 2)

#Transformación en base de datos
data_scale <- as.data.frame(x_norm_scale)

#Sustracción de las sondas
data_scale <- data_scale[grep("A_2|A_3", data_scale$ProbeName) , ]
data_scale <- data_scale[ ,-c(1:3,5)]

#Agregación de la base de datos por método de la media en función de ProbeName
data_agreg_scale <- aggregate(data_scale, by = list(data_scale$ProbeName), FUN = mean, na.rm = TRUE)
data_agreg_scale <- data_agreg_scale[ ,-c(2)]
colnames(data_agreg_scale)[1] = "ID"

#Anotación de la base de datos
data_anot_scale <- data_agreg_scale %>% inner_join(., featuredata, by = "ID")

#Eliminación de las sondas que no tienen correlación con ningún gen
data_anot_scale_noNA <- data_anot_scale[!data_anot_scale$GENE_SYMBOL == "",]

#Segunda agregación de la base de datos por método de la media en función de GENE_SYMBOL
data_final_scale <- aggregate(data_anot_scale_noNA[, -c(1,131)], by = list(data_anot_scale_noNA$GENE_SYMBOL), FUN = mean, na.rm = TRUE)
colnames(data_final_scale)[1] = "Gene"
rownames(data_final_scale) <- NULL
data_final_scale <- column_to_rownames(data_final_scale, var = "Gene")

#Eliminación de las filas que contienen algún valor negativo
data_final_scale[data_final_scale < 0] <- NA
data_final_scale <- na.omit(data_final_scale)

#Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_scale <- prcomp(t(data_final_scale), scale = TRUE)
plot(pca_scale$x[, 1], pca_scale$x[, 2])
pca.var_scale <- pca_scale$sdev^2
pca.var.per_scale <- round(pca.var_scale/sum(pca.var_scale)*100, 1)
barplot(pca.var.per_scale, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_scale <- data.frame(Sample = rownames(pca_scale$x), X = pca_scale$x[, 1], Y = pca_scale$x[, 2])
ggplot(data = pca.data_scale, aes(x = pca_scale$x[, 1], y = pca_scale$x[, 2], label = Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_scale[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_scale[2], "%", sep = "")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_scale <- pca_scale$rotation[, 1]
gene_score_scale <- abs(loading_score_scale)
gene_score_ranked_scale <- sort(gene_score_scale, decreasing = TRUE)
top_10_genes_scale <- names(gene_score_ranked_scale[1:10])
top_10_genes_scale
pca_scale$rotation[top_10_genes_scale, 1]
pca.data2_scale <- data.frame(Sample = rownames(pca_scale$x), X = pca_scale$x[, 1], Y = pca_scale$x[, 2], Condition = diagnostic)
PCA_scale <- ggplot(data = pca.data2_scale, aes(x = pca_scale$x[, 1], y = pca_scale$x[, 2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_scale[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_scale[2], "%", sep="")) + theme_bw() + ggtitle("SCALE")
PCA_scale

#limma
fitlimma_scale <- lmFit(data_final_scale, design)
fitlimma_scale <- contrasts.fit(fitlimma_scale, contrast)
fitlimma_scale <- eBayes(fitlimma_scale)
toptablelimma_scale <- topTable(fitlimma_scale, sort.by = "logFC", number = Inf, adjust.method = "none", coef = 1)
toptablelimma_scale
results_scale <- decideTests(fitlimma_scale)

#Volcano plot
# Si el gen cumple logFC > 0.5 y  P.Value < 0.05 se denomina como "UP" 
toptablelimma_scale$diffexpressed[toptablelimma_scale$logFC > 0.5 & toptablelimma_scale$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_scale$diffexpressed[toptablelimma_scale$logFC < -0.5 & toptablelimma_scale$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_scale$diffexpressed[is.na(toptablelimma_scale$diffexpressed)] <- "NO"
#Función gráfica Volcano
vp_scale <- ggplot(data=toptablelimma_scale, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col = "black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("SCALE") + theme(plot.title = element_text(hjust=0.5))
vp_scale

#mapa de calor
#Selección de las 20 primeros genes que mayor diferencia de expresión presentan
names_scale <- rownames(toptablelimma_scale[1:20,])
df_scale_short <- data_final_scale[c(names_scale),c(1:114)]
#Función del mapa de calor
heatmap_scale <- pheatmap(df_scale_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Scale", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row")



#Método de normalización "cyclicloess"

x_norm_cyclicloess <- normalizeBetweenArrays(x_fondo, method = "cyclicloess")

#Control distribución después de normalizar y eliminar fondo 
plot_norm_cyclicloess <- plotDensities(x_norm_cyclicloess, col = "black", legend = FALSE)
plotMD1_norm_cyclicloess <- plotMD(x_norm_cyclicloess, column = 1) + abline(h = 0, col = "red", lty = 2, lwd = 2)
plotMD2_norm_cyclicloess <- plotMD(x_norm_cyclicloess, column = 2) + abline(h = 0, col = "red", lty = 2, lwd = 2)

#Transformación en base de datos
data_cyclicloess <- as.data.frame(x_norm_cyclicloess)

#Sustracción de las sondas
data_cyclicloess <- data_cyclicloess[grep("A_2|A_3", data_cyclicloess$ProbeName) , ]
data_cyclicloess <- data_cyclicloess[ ,-c(1:3, 5)]

#Agregación de la base de datos por método de la media en función de ProbeName
data_agreg_cyclicloess <- aggregate(data_cyclicloess, by = list(data_cyclicloess$ProbeName), FUN = mean, na.rm = TRUE)
data_agreg_cyclicloess <- data_agreg_cyclicloess[ ,-c(2)]
colnames(data_agreg_cyclicloess)[1] = "ID"

#Anotación de la base de datos
data_anot_cyclicloess <- data_agreg_cyclicloess %>% inner_join(., featuredata, by = "ID")

#Eliminación de las sondas que no tienen correlación con ningún gen
data_anot_cyclicloess_noNA <- data_anot_cyclicloess[!data_anot_cyclicloess$GENE_SYMBOL == "",]

#Segunda agregación de la base de datos por método de la media en función de GENE_SYMBOL
data_final_cyclicloess <- aggregate(data_anot_cyclicloess_noNA[, -c(1, 131)], by = list(data_anot_cyclicloess_noNA$GENE_SYMBOL), FUN = mean, na.rm = TRUE)
colnames(data_final_cyclicloess)[1] = "Gene"
rownames(data_final_cyclicloess) <- NULL
data_final_cyclicloess <- column_to_rownames(data_final_cyclicloess, var = "Gene")

#Eliminación de las filas que contienen algún valor negativo
data_final_cyclicloess[data_final_cyclicloess < 0] <- NA
data_final_cyclicloess <- na.omit(data_final_cyclicloess)

#Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_cyclicloess <- prcomp(t(data_final_cyclicloess), scale = TRUE)
plot(pca_cyclicloess$x[, 1], pca_cyclicloess$x[, 2])
pca.var_cyclicloess <- pca_cyclicloess$sdev^2
pca.var.per_cyclicloess <- round(pca.var_cyclicloess/sum(pca.var_cyclicloess)*100, 1)
barplot(pca.var.per_cyclicloess, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_cyclicloess <- data.frame(Sample = rownames(pca_cyclicloess$x), X = pca_cyclicloess$x[, 1], Y = pca_cyclicloess$x[, 2])
ggplot(data = pca.data_cyclicloess, aes(x = pca_cyclicloess$x[, 1], y = pca_cyclicloess$x[, 2], label = Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_cyclicloess[1], "%", sep = "")) + ylab(paste("PC2 - ", pca.var.per_cyclicloess[2], "%", sep = "")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_cyclicloess <- pca_cyclicloess$rotation[, 1]
gene_score_cyclicloess <- abs(loading_score_cyclicloess)
gene_score_ranked_cyclicloess <- sort(gene_score_cyclicloess, decreasing=TRUE)
top_10_genes_cyclicloess <- names(gene_score_ranked_cyclicloess[1:10])
top_10_genes_cyclicloess
#para mostrar las puntuaciones
pca_cyclicloess$rotation[top_10_genes_cyclicloess, 1]
#color grupos
#f <- factor(metadatos$diagnosis)
pca.data2_cyclicloess <- data.frame(Sample = rownames(pca_cyclicloess$x), X = pca_cyclicloess$x[,1], Y = pca_cyclicloess$x[,2], Condition = diagnostic)
PCA_cyclicloess <- ggplot(data = pca.data2_cyclicloess, aes(x = pca_cyclicloess$x[, 1], y = pca_cyclicloess$x[, 2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_cyclicloess[1], "%", sep = "")) + ylab(paste("PC2 - ", pca.var.per_cyclicloess[2], "%", sep = "")) + theme_bw() + ggtitle("CYCLICLOESS")
PCA_cyclicloess

#limma
fitlimma_cyclicloess <- lmFit(data_final_cyclicloess, design)
fitlimma_cyclicloess <- contrasts.fit(fitlimma_cyclicloess, contrast)
fitlimma_cyclicloess <- eBayes(fitlimma_cyclicloess)
toptablelimma_cyclicloess <- topTable(fitlimma_cyclicloess, sort.by = "logFC", number = Inf, adjust.method = "none", coef = 1)
toptablelimma_cyclicloess
results_cyclicloess <- decideTests(fitlimma_cyclicloess)

#Volcano plot
# Si el gen cumple logFC > 0.5 y  P.Value < 0.05 se denomina como "UP" 
toptablelimma_cyclicloess$diffexpressed[toptablelimma_cyclicloess$logFC > 0.5 & toptablelimma_cyclicloess$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_cyclicloess$diffexpressed[toptablelimma_cyclicloess$logFC < -0.5 & toptablelimma_cyclicloess$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_cyclicloess$diffexpressed[is.na(toptablelimma_cyclicloess$diffexpressed)] <- "NO"
#Función gráfica Volcano
vp_cyclicloess <- ggplot(data=toptablelimma_cyclicloess, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col="black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(hjust = 0.5))
vp_cyclicloess

#mapa de calor
#Selección de las 20 primeros genes que mayor diferencia de expresión presentan
names_cyclicloess <- rownames(toptablelimma_cyclicloess[1:20,])
df_cyclicloess_short <- data_final_cyclicloess[c(names_cyclicloess), c(1:114)]
#Función del mapa de calor
heatmap_cyclicloess <- pheatmap(df_cyclicloess_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Cyclicloess", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row")


#Diagrama de Venn
UP_none <- row.names(toptablelimma_none)[which(toptablelimma_none$diffexpressed == "UP")]
UP_quantile <- row.names(toptablelimma_quantile)[which(toptablelimma_quantile$diffexpressed == "UP")]
UP_scale <- row.names(toptablelimma_scale)[which(toptablelimma_scale$diffexpressed == "UP")]
UP_cyclicloess <- row.names(toptablelimma_cyclicloess)[which(toptablelimma_cyclicloess$diffexpressed == "UP")]

DOWN_none <- row.names(toptablelimma_none)[which(toptablelimma_none$diffexpressed == "DOWN")]
DOWN_quantile <- row.names(toptablelimma_quantile)[which(toptablelimma_quantile$diffexpressed == "DOWN")]
DOWN_scale <- row.names(toptablelimma_scale)[which(toptablelimma_scale$diffexpressed == "DOWN")]
DOWN_cyclicloess <- row.names(toptablelimma_cyclicloess)[which(toptablelimma_cyclicloess$diffexpressed == "DOWN")]

length(UP_none)
length(UP_quantile)
length(UP_scale)
length(UP_cyclicloess)

length(DOWN_none)
length(DOWN_quantile)
length(DOWN_scale)
length(DOWN_cyclicloess)

lista_UP <- list(UP_none, UP_quantile, UP_scale, UP_cyclicloess)
lista_DOWN <- list(DOWN_none, DOWN_quantile, DOWN_scale, DOWN_cyclicloess)

ggVennDiagram(lista_UP, category.names = c("NONE","QUANTILE","SCALE", "CYCLICLOESS"), label = "count", set.size = 2, label_alpha = 0) + scale_fill_distiller(palette = "RdBu") + scale_x_continuous(expand = expansion(mult = .2)) + ggtitle("SOBREEXPRESADOS") + theme(plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
ggVennDiagram(lista_DOWN, category.names = c("NONE","QUANTILE","SCALE", "CYCLICLOESS"), label = "count", set.size = 2, label_alpha = 0) + scale_fill_distiller(palette = "RdBu") + scale_x_continuous(expand = expansion(mult = .2)) + ggtitle("INFRAEXPRESADOS") + theme(plot.title = element_text(face = "bold", size = 30, hjust = 0.5))

UP10_none <- UP_none[1:10]
UP100_none <- UP_none[1:100]
DOWN10_none <- DOWN_none[1:10]
DOWN100_none <- DOWN_none[1:100]

UP10_quantile <- UP_quantile[1:10]
UP100_quantile <- UP_quantile[1:100]
DOWN10_quantile <- DOWN_quantile[1:10]
DOWN100_quantile <- DOWN_quantile[1:100]

UP10_scale <- UP_scale[1:10]
UP100_scale <- UP_scale[1:100]
DOWN10_scale <- DOWN_scale[1:10]
DOWN100_scale <- DOWN_scale[1:100]

UP10_cyclicloess <- UP_cyclicloess[1:10]
UP100_cyclicloess <- UP_cyclicloess[1:100]
DOWN10_cyclicloess <- DOWN_cyclicloess[1:10]
DOWN100_cyclicloess <- DOWN_cyclicloess[1:100]

lista_UP10 <- list(UP10_none, UP10_quantile, UP10_scale, UP10_cyclicloess)
lista_DOWN10 <- list(DOWN10_none, DOWN10_quantile, DOWN10_scale, DOWN10_cyclicloess)

ggVennDiagram(lista_UP10, category.names = c("NONE","QUANTILE","SCALE", "CYCLICLOESS"), label = "count", set.size = 2, label_alpha = 0) + scale_fill_distiller(palette = "RdBu") + scale_x_continuous(expand = expansion(mult = .2)) + ggtitle("10 SOBREEXPRESADOS") + theme(plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
ggVennDiagram(lista_DOWN10, category.names = c("NONE","QUANTILE","SCALE", "CYCLICLOESS"), label = "count", set.size = 2, label_alpha = 0) + scale_fill_distiller(palette = "RdBu") + scale_x_continuous(expand = expansion(mult = .2)) + ggtitle("10 INFRAEXPRESADOS") + theme(plot.title = element_text(face = "bold", size = 30, hjust = 0.5))

lista_UP100 <- list(UP100_none, UP100_quantile, UP100_scale, UP100_cyclicloess)
lista_DOWN100 <- list(DOWN100_none, DOWN100_quantile, DOWN100_scale, DOWN100_cyclicloess)

ggVennDiagram(lista_UP100, category.names = c("NONE","QUANTILE","SCALE", "CYCLICLOESS"), label = "count", set.size = 2, label_alpha = 0) + scale_fill_distiller(palette = "RdBu") + scale_x_continuous(expand = expansion(mult = .2)) + ggtitle("100 SOBREEXPRESADOS") + theme(plot.title = element_text(face = "bold", size = 30, hjust = 0.5))
ggVennDiagram(lista_DOWN100, category.names = c("NONE","QUANTILE","SCALE", "CYCLICLOESS"), label = "count", set.size = 2, label_alpha = 0) + scale_fill_distiller(palette = "RdBu") + scale_x_continuous(expand = expansion(mult = .2)) + ggtitle("100 INFRAEXPRESADOS") + theme(plot.title = element_text(face = "bold", size = 30, hjust = 0.5))


#Figura PCA
PCA_none2 <- ggplot(data = pca.data_none2, aes(x = pca_none$x[,1], y = pca_none$x[,2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep = "")) + theme_bw() + ggtitle("NONE") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")
PCA_quantile2 <- ggplot(data = pca.data2_quantile, aes(x = pca_quantile$x[,1], y = pca_quantile$x[,2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_quantile[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_quantile[2], "%", sep = "")) + theme_bw() + ggtitle("QUANTILE") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")
PCA_scale2 <- ggplot(data = pca.data2_scale, aes(x = pca_scale$x[,1], y = pca_scale$x[,2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_scale[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_scale[2], "%", sep = "")) + theme_bw() + ggtitle("SCALE") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")
PCA_cyclicloess2 <- ggplot(data = pca.data2_cyclicloess, aes(x = pca_cyclicloess$x[,1], y = pca_cyclicloess$x[,2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_cyclicloess[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_cyclicloess[2], "%", sep = "")) + theme_bw() + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")

ext_leyenda_PCA <- ggplot(data = pca.data_none2, aes(x = pca_none$x[,1], y = pca_none$x[,2], label = Sample, color = Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep = "")) + theme_bw() + ggtitle("NONE") + theme(legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12))
ext_leyenda_PCA <- ext_leyenda_PCA
leyenda_PCA <- get_legend(ext_leyenda_PCA)

PCA_grupo <- plot_grid(PCA_none2, NULL, PCA_quantile2, PCA_scale2, NULL, PCA_cyclicloess2, ncol = 3, nrow = 2, labels = NULL, rel_widths = c(1, 0.06, 1))
PCA_figura <- plot_grid(leyenda_PCA, PCA_grupo, ncol=1, rel_heights = c(1, 8))
PCA_figura


#Figuras mapa de calor
heatmap_none2 <- pheatmap(df_none_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "NONE", fontsize = 14, colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)
heatmap_quantile2 <- pheatmap(df_quantile_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE,show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "QUANTILE", fontsize = 14, colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)
heatmap_scale2 <- pheatmap(df_scale_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE,show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "SCALE", fontsize = 14, colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)
heatmap_cyclicloess2 <- pheatmap(df_cyclicloess_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE,show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "CYCLICLOESS", fontsize = 14, colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)


#Figura volcano
vp_none2 <- ggplot(data = toptablelimma_none, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col="black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("NONE") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")
vp_quantile2 <- ggplot(data = toptablelimma_quantile, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col="black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("QUANTILE") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")
vp_scale2 <- ggplot(data = toptablelimma_scale, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col="black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("SCALE") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")
vp_cyclicloess2 <- ggplot(data = toptablelimma_cyclicloess, aes(x = logFC, y = -log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col="black") + geom_hline(yintercept = -log10(0.05), col = "black") + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + theme(legend.position = "none")

ext_leyenda_vp <- ggplot(data = toptablelimma_cyclicloess, aes(x = logFC, y = -log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values = c("red", "black", "blue")) + geom_vline(xintercept = c(-0.5, 0.5), col="black") + geom_hline(yintercept = -log10(0.05), col="black") + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12))
leyenda_vp <- get_legend(ext_leyenda_vp)

volcano_grupo <- plot_grid(vp_none2, NULL, vp_quantile2, vp_scale2, NULL, vp_cyclicloess2, ncol = 3, nrow = 2, labels = NULL, rel_widths = c(1, 0.06, 1))
volcanos_figura <- plot_grid(leyenda_vp, volcano_grupo, ncol=1, rel_heights = c(1, 8))
volcanos_figura

