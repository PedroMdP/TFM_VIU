#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Adding_annotation
#https://rnnh.github.io/bioinfo-notebook/docs/DE_analysis_edgeR_script.html

#How do you analyze microarray data?
#1. feature extraction.
#2. quality control.
#3. normalisation.
#4. differential expression analysis.
#5. biological interpretation of the results.
#6. submission of data to a public database.


library(Biobase)
library(limma)
library(VennDiagram)
library(ggplot2)
library(gplots)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)


#Establecimiento del directorio de trabajo
setwd("/Users/pedro/Documents/VIU/0 - Trabajo Fin de Máster/Experimentación/Lecturas prueba")

#Lectura metadatos, creación de grupos y tabla "diseño"
metadata <- read.delim("metadata.csv", sep=";", row.name = 1)
group <- paste(metadata$diagnosis, sep=";")
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
featuredata <- read.table("GPL13497-9755.txt", sep="\t", header = TRUE, fill = TRUE, quote="\"")
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
background <- boxplot(data.frame(log2(x$Eb)), main="Background")
plotMD1_no_norm <- plotMD(x, column = 1) + abline(h=0, col="red", lty=2, lwd=2)
plotMD2_no_norm <- plotMD(x, column = 2) + abline(h=0, col="red", lty=2, lwd=2)
#Intentar que funcione para añadir a los anteriores gráficos
#boxplot(exprs(e.raw), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

#set parameters and draw the plot
#dev.new(width=4+dim(e.raw)[[2]]/5, height=6)
#par(mar=c(2+round(max(nchar(sampleNames(e.raw)))/2),4,2,1))
#title <- paste ("GSE131761", '/', annotation(e.raw), " selected samples", sep ='')
#boxplot(exprs(e.raw), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
#plotMD(e.raw, column = 1)
#plotMD(e.raw, column = 2)

#Corrección de fondo y normalización
x_fondo <- backgroundCorrect(x, method="normexp")


#Método de normalización "none"

x_norm_none <- normalizeBetweenArrays(x_fondo, method="none")

#Control distribución después de la eliminación de fondo y de la normalización 
plot_norm_none <- plotDensities(x_norm_none, col = "black", legend = FALSE)
background_norm_none <- boxplot(data.frame(log2(x_norm_none$E)),main="Background")
plotMD1_norm_none <- plotMD(x_norm_none, column = 1) + abline(h=0, col="red", lty=2, lwd=2)
plotMD2_norm_none <- plotMD(x_norm_none, column = 2) + abline(h=0, col="red", lty=2, lwd=2)

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
data_anot_none_noNA <- data_anot_none[!data_anot_none$GENE_SYMBOL=="",]

#Segunda agregación de la base de datos por método de la media en función de GENE_SYMBOL
data_final_none <- aggregate(data_anot_none_noNA[, -c(1,131)], by = list(data_anot_none_noNA$GENE_SYMBOL), FUN = mean, na.rm = TRUE)
colnames(data_final_none)[1] = "Gene"
rownames(data_final_none) <- NULL
data_final_none <- column_to_rownames(data_final_none, var = "Gene")

#Eliminación de las filas que contienen algún valor negativo (Comprobar) => añadir remove.zeros = FALSE en DGEList
data_final_none[data_final_none < 0] <- NA
data_final_none <- na.omit(data_final_none)

#(https://www.youtube.com/watch?v=0Jp4gsfOLMs) Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_none <- prcomp(t(data_final_none), scale=TRUE)
plot(pca_none$x[,1], pca_none$x[,2])
pca.var_none <- pca_none$sdev^2
pca.var.per_none <- round(pca.var_none/sum(pca.var_none)*100, 1)
barplot(pca.var.per_none, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_none <- data.frame(Sample=rownames(pca_none$x), X=pca_none$x[,1], Y=pca_none$x[,2])
ggplot(data=pca.data_none, aes(x=pca_none$x[,1], y=pca_none$x[,2], label=Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep="")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_none <- pca_none$rotation[,1]
gene_score_none <- abs(loading_score_none)
gene_score_ranked_none <- sort(gene_score_none, decreasing=TRUE)
top_10_genes_none <- names(gene_score_ranked_none[1:10])
top_10_genes_none
#para mostrar las puntuaciones
pca_none$rotation[top_10_genes_none,1]
#color grupos
#f <- factor(metadatos$diagnosis)
pca.data_none2 <- data.frame(Sample=rownames(pca_none$x), X=pca_none$x[,1], Y=pca_none$x[,2], Condition=diagnostic)
PCA_none <- ggplot(data=pca.data_none2, aes(x=pca_none$x[,1], y=pca_none$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep="")) + theme_bw() + ggtitle("NONE")


#limma
fitlimma_none <- lmFit(data_final_none, design)
fitlimma_none <- contrasts.fit(fitlimma_none, contrast)
fitlimma_none <- eBayes(fitlimma_none)
toptablelimma_none <- topTable(fitlimma_none, sort.by = "logFC", number=Inf, adjust.method = "none", coef = 1)
toptablelimma_none
results_none <- decideTests(fitlimma_none)
#vennDiagram(results_none)
length(which(toptablelimma_none$adj.P.Val < 0.05))
siglimma_none <- subset(toptablelimma_none, adj.P.Val<0.05)
siglimma_up_none <- rownames(subset(siglimma_none,logFC>0))
siglimma_dn_none <- rownames(subset(siglimma_none,logFC<0))
length(siglimma_up_none)
length(siglimma_dn_none)

#Volcano plot
# Si el gen cumple logFC > 0.5 y P.Value < 0.05 se denomina como "UP" 
toptablelimma_none$diffexpressed[toptablelimma_none$logFC > 0.5 & toptablelimma_none$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_none$diffexpressed[toptablelimma_none$logFC < -0.5 & toptablelimma_none$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_none$diffexpressed[is.na(toptablelimma_none$diffexpressed)] <- "NO_none"
#Función gráfica Volcano
vp_none <- ggplot(data=toptablelimma_none, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("NONE") + theme(plot.title = element_text(hjust=0.5))
vp_none

#Mapa de calor
#Selección de los 20 primeros genes que mayor diferencia de expresión presentan
names_none <- rownames(toptablelimma_none[1:20,])
df_none_short <- data_final_none[c(names_none),c(1:114)]
#Anotación de las muestras
annotation_df <- metadata[1:114,4,drop = FALSE]
#Función mapa de calor
heatmap_none <- pheatmap(df_none_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "None", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row")


#Método de normalización "quantile"

x_norm_quantile <- normalizeBetweenArrays(x_fondo, method="quantile")

#Control distribución después de normalizar y eliminar fondo 
plot_norm_quantile <- plotDensities(x_norm_quantile, col = "black", legend = FALSE)
plotMD1_norm_quantile <- plotMD(x_norm_quantile, column = 1) + abline(h=0, col="red", lty=2, lwd=2)
plotMD2_norm_quantile <- plotMD(x_norm_quantile, column = 2) + abline(h=0, col="red", lty=2, lwd=2)

#Transformación en base de datos
data_quantile <- as.data.frame(x_norm_quantile)

#Sustracción de las sondas
data_quantile <- data_quantile[grep("A_2|A_3", data_quantile$ProbeName) , ]
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

#Eliminación de las filas que contienen algún valor negativo (Comprobar) => añadir remove.zeros = FALSE en DGEList
data_final_quantile[data_final_quantile < 0] <- NA
data_final_quantile <- na.omit(data_final_quantile)

#(https://www.youtube.com/watch?v=0Jp4gsfOLMs) Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_quantile <- prcomp(t(data_final_quantile), scale=TRUE)
plot(pca_quantile$x[,1], pca_quantile$x[,2])
pca.var_quantile <- pca_quantile$sdev^2
pca.var.per_quantile <- round(pca.var_quantile/sum(pca.var_quantile)*100, 1)
barplot(pca.var.per_quantile, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_quantile <- data.frame(Sample=rownames(pca_quantile$x), X=pca_quantile$x[,1], Y=pca_quantile$x[,2])
ggplot(data=pca.data_quantile, aes(x=pca_quantile$x[,1], y=pca_quantile$x[,2], label=Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_quantile[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_quantile[2], "%", sep="")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_quantile <- pca_quantile$rotation[,1]
gene_score_quantile <- abs(loading_score_quantile)
gene_score_ranked_quantile <- sort(gene_score_quantile, decreasing=TRUE)
top_10_genes_quantile <- names(gene_score_ranked_quantile[1:10])
top_10_genes_quantile
#para mostrar las puntuaciones
pca_quantile$rotation[top_10_genes_quantile,1]
#color grupos
#f <- factor(metadatos$diagnosis)
pca.data2_quantile <- data.frame(Sample=rownames(pca_quantile$x), X=pca_quantile$x[,1], Y=pca_quantile$x[,2], Condition=diagnostic)
PCA_quantile <- ggplot(data=pca.data2_quantile, aes(x=pca_quantile$x[,1], y=pca_quantile$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_quantile[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_quantile[2], "%", sep="")) + theme_bw() + ggtitle("QUANTILE")

#limma
fitlimma_quantile <- lmFit(data_final_quantile, design)
fitlimma_quantile <- contrasts.fit(fitlimma_quantile, contrast)
fitlimma_quantile <- eBayes(fitlimma_quantile)
toptablelimma_quantile <- topTable(fitlimma_quantile, sort.by = "logFC", number=Inf, adjust.method = "none", coef = 1)
toptablelimma_quantile
results_quantile <- decideTests(fitlimma_quantile)
#vennDiagram(results_quantile)
length(which(toptablelimma_quantile$adj.P.Val < 0.05))
siglimma_quantile <- subset(toptablelimma_quantile, adj.P.Val<0.05)
siglimma_up_quantile <- rownames(subset(siglimma_quantile,logFC>0))
siglimma_dn_quantile <- rownames(subset(siglimma_quantile,logFC<0))
length(siglimma_up_quantile)
length(siglimma_dn_quantile)

#Volcano plot
# Si el gen cumple logFC > 0.5 y  P.Value < 0.05 se denomina como "UP" 
toptablelimma_quantile$diffexpressed[toptablelimma_quantile$logFC > 0.5 & toptablelimma_quantile$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_quantile$diffexpressed[toptablelimma_quantile$logFC < -0.5 & toptablelimma_quantile$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_quantile$diffexpressed[is.na(toptablelimma_quantile$diffexpressed)] <- "NO_none"
#Función gráfica Volcano
vp_quantile <- ggplot(data=toptablelimma_quantile, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("QUANTILE") + theme(plot.title = element_text(hjust=0.5))
vp_quantile

#mapa de calor
#Selección de las 20 primeros genes que mayor diferencia de expresión presentan
names_quantile <- rownames(toptablelimma_quantile[1:20,])
df_quantile_short <- data_final_quantile[c(names_quantile),c(1:114)]
#Función del mapa de calor
heatmap_quantile <- pheatmap(df_quantile_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Quantile", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row")



#MÉTODO NORMALIZACIÓN "SCALE"

x_norm_scale <- normalizeBetweenArrays(x_fondo, method="scale")

#Control distribución después de normalizar y eliminar fondo 
plot_norm_scale <- plotDensities(x_norm_scale, col = "black", legend = FALSE)
plotMD1_norm_scale <- plotMD(x_norm_scale, column = 1) + abline(h=0, col="red", lty=2, lwd=2)
plotMD2_norm_scale <- plotMD(x_norm_scale, column = 2) + abline(h=0, col="red", lty=2, lwd=2)

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
data_anot_scale_noNA <- data_anot_scale[!data_anot_scale$GENE_SYMBOL=="",]

#Segunda agregación de la base de datos por método de la media en función de GENE_SYMBOL
data_final_scale <- aggregate(data_anot_scale_noNA[, -c(1,131)], by = list(data_anot_scale_noNA$GENE_SYMBOL), FUN = mean, na.rm = TRUE)
colnames(data_final_scale)[1] = "Gene"
rownames(data_final_scale) <- NULL
data_final_scale <- column_to_rownames(data_final_scale, var = "Gene")

#Eliminación de las filas que contienen algún valor negativo (Comprobar) => añadir remove.zeros = FALSE en DGEList
data_final_scale[data_final_scale < 0] <- NA
data_final_scale <- na.omit(data_final_scale)

#(https://www.youtube.com/watch?v=0Jp4gsfOLMs) Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_scale <- prcomp(t(data_final_scale), scale=TRUE)
plot(pca_scale$x[,1], pca_scale$x[,2])
pca.var_scale <- pca_scale$sdev^2
pca.var.per_scale <- round(pca.var_scale/sum(pca.var_scale)*100, 1)
barplot(pca.var.per_scale, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_scale <- data.frame(Sample=rownames(pca_scale$x), X=pca_scale$x[,1], Y=pca_scale$x[,2])
ggplot(data=pca.data_scale, aes(x=pca_scale$x[,1], y=pca_scale$x[,2], label=Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_scale[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_scale[2], "%", sep="")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_scale <- pca_scale$rotation[,1]
gene_score_scale <- abs(loading_score_scale)
gene_score_ranked_scale <- sort(gene_score_scale, decreasing=TRUE)
top_10_genes_scale <- names(gene_score_ranked_scale[1:10])
top_10_genes_scale
#para mostrar las puntuaciones
pca_scale$rotation[top_10_genes_scale,1]
#color grupos
#f <- factor(metadatos$diagnosis)
pca.data2_scale <- data.frame(Sample=rownames(pca_scale$x), X=pca_scale$x[,1], Y=pca_scale$x[,2], Condition=diagnostic)
PCA_scale <- ggplot(data=pca.data2_scale, aes(x=pca_scale$x[,1], y=pca_scale$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_scale[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_scale[2], "%", sep="")) + theme_bw() + ggtitle("SCALE")

#limma
fitlimma_scale <- lmFit(data_final_scale, design)
fitlimma_scale <- contrasts.fit(fitlimma_scale, contrast)
fitlimma_scale <- eBayes(fitlimma_scale)
toptablelimma_scale <- topTable(fitlimma_scale, sort.by = "logFC", number=Inf, adjust.method = "none", coef = 1)
toptablelimma_scale
results_scale <- decideTests(fitlimma_scale)
#vennDiagram(results_scale)
length(which(toptablelimma_scale$adj.P.Val < 0.05))
siglimma_scale <- subset(toptablelimma_scale, adj.P.Val<0.05)
siglimma_up_scale <- rownames(subset(siglimma_scale,logFC>0))
siglimma_dn_scale <- rownames(subset(siglimma_scale,logFC<0))
length(siglimma_up_scale)
length(siglimma_dn_scale)

#Volcano plot
# Si el gen cumple logFC > 0.5 y  P.Value < 0.05 se denomina como "UP" 
toptablelimma_scale$diffexpressed[toptablelimma_scale$logFC > 0.5 & toptablelimma_scale$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_scale$diffexpressed[toptablelimma_scale$logFC < -0.5 & toptablelimma_scale$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_scale$diffexpressed[is.na(toptablelimma_scale$diffexpressed)] <- "NO"
#Función gráfica Volcano
vp_scale <- ggplot(data=toptablelimma_scale, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("SCALE") + theme(plot.title = element_text(hjust=0.5))
vp_scale

#mapa de calor
#Selección de las 20 primeros genes que mayor diferencia de expresión presentan
names_scale <- rownames(toptablelimma_scale[1:20,])
df_scale_short <- data_final_scale[c(names_scale),c(1:114)]
#Función del mapa de calor
heatmap_scale <- pheatmap(df_scale_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Scale", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row")



#Método de normalización "cyclicloess"

x_norm_cyclicloess <- normalizeBetweenArrays(x_fondo, method="cyclicloess")

#Control distribución después de normalizar y eliminar fondo 
plot_norm_cyclicloess <- plotDensities(x_norm_cyclicloess, col = "black", legend = FALSE)
plotMD1_norm_cyclicloess <- plotMD(x_norm_cyclicloess, column = 1) + abline(h=0, col="red", lty=2, lwd=2)
plotMD2_norm_cyclicloess <- plotMD(x_norm_cyclicloess, column = 2) + abline(h=0, col="red", lty=2, lwd=2)

#Transformación en base de datos
data_cyclicloess <- as.data.frame(x_norm_cyclicloess)

#Sustracción de las sondas
data_cyclicloess <- data_cyclicloess[grep("A_2|A_3", data_cyclicloess$ProbeName) , ]
data_cyclicloess <- data_cyclicloess[ ,-c(1:3,5)]

#Agregación de la base de datos por método de la media en función de ProbeName
data_agreg_cyclicloess <- aggregate(data_cyclicloess, by = list(data_cyclicloess$ProbeName), FUN = mean, na.rm = TRUE)
data_agreg_cyclicloess <- data_agreg_cyclicloess[ ,-c(2)]
colnames(data_agreg_cyclicloess)[1] = "ID"

#Anotación de la base de datos
data_anot_cyclicloess <- data_agreg_cyclicloess %>% inner_join(., featuredata, by = "ID")

#Eliminación de las sondas que no tienen correlación con ningún gen
data_anot_cyclicloess_noNA <- data_anot_cyclicloess[!data_anot_cyclicloess$GENE_SYMBOL=="",]

#Segunda agregación de la base de datos por método de la media en función de GENE_SYMBOL
data_final_cyclicloess <- aggregate(data_anot_cyclicloess_noNA[, -c(1,131)], by = list(data_anot_cyclicloess_noNA$GENE_SYMBOL), FUN = mean, na.rm = TRUE)
colnames(data_final_cyclicloess)[1] = "Gene"
rownames(data_final_cyclicloess) <- NULL
data_final_cyclicloess <- column_to_rownames(data_final_cyclicloess, var = "Gene")

#Eliminación de las filas que contienen algún valor negativo (Comprobar) => añadir remove.zeros = FALSE en DGEList
data_final_cyclicloess[data_final_cyclicloess < 0] <- NA
data_final_cyclicloess <- na.omit(data_final_cyclicloess)

#(https://www.youtube.com/watch?v=0Jp4gsfOLMs) Cálculo de PCA -> representación -> cálculo de la variación -> cálculo del porcentaje de variación -> representación del porcentaje de variación
pca_cyclicloess <- prcomp(t(data_final_cyclicloess), scale=TRUE)
plot(pca_cyclicloess$x[,1], pca_cyclicloess$x[,2])
pca.var_cyclicloess <- pca_cyclicloess$sdev^2
pca.var.per_cyclicloess <- round(pca.var_cyclicloess/sum(pca.var_cyclicloess)*100, 1)
barplot(pca.var.per_cyclicloess, main = "Gráfica", xlab = "Componente Principal", ylab = "Porcentaje variación")

#PCA con más información
pca.data_cyclicloess <- data.frame(Sample=rownames(pca_cyclicloess$x), X=pca_cyclicloess$x[,1], Y=pca_cyclicloess$x[,2])
ggplot(data=pca.data_cyclicloess, aes(x=pca_cyclicloess$x[,1], y=pca_cyclicloess$x[,2], label=Sample)) + geom_text() + xlab(paste("PC1 - ", pca.var.per_cyclicloess[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_cyclicloess[2], "%", sep="")) + theme_bw() + ggtitle("Gráfica PCA")
loading_score_cyclicloess <- pca_cyclicloess$rotation[,1]
gene_score_cyclicloess <- abs(loading_score_cyclicloess)
gene_score_ranked_cyclicloess <- sort(gene_score_cyclicloess, decreasing=TRUE)
top_10_genes_cyclicloess <- names(gene_score_ranked_cyclicloess[1:10])
top_10_genes_cyclicloess
#para mostrar las puntuaciones
pca_cyclicloess$rotation[top_10_genes_cyclicloess,1]
#color grupos
#f <- factor(metadatos$diagnosis)
pca.data2_cyclicloess <- data.frame(Sample=rownames(pca_cyclicloess$x), X=pca_cyclicloess$x[,1], Y=pca_cyclicloess$x[,2], Condition=diagnostic)
PCA_cyclicloess <- ggplot(data=pca.data2_cyclicloess, aes(x=pca_cyclicloess$x[,1], y=pca_cyclicloess$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_cyclicloess[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_cyclicloess[2], "%", sep="")) + theme_bw() + ggtitle("CYCLICLOESS")

#limma
fitlimma_cyclicloess <- lmFit(data_final_cyclicloess, design)
fitlimma_cyclicloess <- contrasts.fit(fitlimma_cyclicloess, contrast)
fitlimma_cyclicloess <- eBayes(fitlimma_cyclicloess)
toptablelimma_cyclicloess <- topTable(fitlimma_cyclicloess, sort.by = "logFC", number=Inf, adjust.method = "none", coef = 1)
toptablelimma_cyclicloess
results_cyclicloess <- decideTests(fitlimma_cyclicloess)
#vennDiagram(results_cyclicloess)
length(which(toptablelimma_cyclicloess$adj.P.Val < 0.05))
siglimma_cyclicloess <- subset(toptablelimma_cyclicloess, adj.P.Val<0.05)
siglimma_up_cyclicloess <- rownames(subset(siglimma_cyclicloess,logFC>0))
siglimma_dn_cyclicloess <- rownames(subset(siglimma_cyclicloess,logFC<0))
length(siglimma_up_cyclicloess)
length(siglimma_dn_cyclicloess)

#Volcano plot
# Si el gen cumple logFC > 0.5 y  P.Value < 0.05 se denomina como "UP" 
toptablelimma_cyclicloess$diffexpressed[toptablelimma_cyclicloess$logFC > 0.5 & toptablelimma_cyclicloess$P.Value < 0.05] <- "UP"
# Si el gen cumple logFC < -0.5 y P.Value < 0.05 se denomina como "DOWN"
toptablelimma_cyclicloess$diffexpressed[toptablelimma_cyclicloess$logFC < -0.5 & toptablelimma_cyclicloess$P.Value < 0.05] <- "DOWN"
# Si el gen cumple logFC > 0.5 o < -0.5 and P.Value < 0.05 se denomina como "NO"
toptablelimma_cyclicloess$diffexpressed[is.na(toptablelimma_cyclicloess$diffexpressed)] <- "NO"
#Función gráfica Volcano
vp_cyclicloess <- ggplot(data=toptablelimma_cyclicloess, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(hjust=0.5))
vp_cyclicloess

#mapa de calor
#Selección de las 20 primeros genes que mayor diferencia de expresión presentan
names_cyclicloess <- rownames(toptablelimma_cyclicloess[1:20,])
df_cyclicloess_short <- data_final_cyclicloess[c(names_cyclicloess),c(1:114)]
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

venn(list("NONE"=UP_none, "QUANTILE"=UP_quantile, "SCALE"=UP_scale, "CYCLICLOESS"=UP_cyclicloess))
venn(list("NONE"=DOWN_none, "QUANTILE"=DOWN_quantile, "SCALE"=DOWN_scale, "CYCLICLOESS"=DOWN_cyclicloess))


UP10_none <- UP_none[1:10]
UP100_none <- UP_none[1:100]
UP1000_none <- UP_none[1:1000]
DOWN10_none <- DOWN_none[1:10]
DOWN100_none <- DOWN_none[1:100]
DOWN1000_none <- DOWN_none[1:1000]

UP10_quantile <- UP_quantile[1:10]
UP100_quantile <- UP_quantile[1:100]
UP1000_quantile <- UP_quantile[1:1000]
DOWN10_quantile <- DOWN_quantile[1:10]
DOWN100_quantile <- DOWN_quantile[1:100]
DOWN1000_quantile <- DOWN_quantile[1:1000]

UP10_scale <- UP_scale[1:10]
UP100_scale <- UP_scale[1:100]
UP1000_scale <- UP_scale[1:1000]
DOWN10_scale <- DOWN_scale[1:10]
DOWN100_scale <- DOWN_scale[1:100]
DOWN1000_scale <- DOWN_scale[1:1000]

UP10_cyclicloess <- UP_cyclicloess[1:10]
UP100_cyclicloess <- UP_cyclicloess[1:100]
UP1000_cyclicloess <- UP_cyclicloess[1:1000]
DOWN10_cyclicloess <- DOWN_cyclicloess[1:10]
DOWN100_cyclicloess <- DOWN_cyclicloess[1:100]
DOWN1000_cyclicloess <- DOWN_cyclicloess[1:1000]

venn(list("NONE"=UP10_none, "QUANTILE"=UP10_quantile, "SCALE"=UP10_scale, "CYCLICLOESS"=UP10_cyclicloess))
venn(list("NONE"=DOWN10_none, "QUANTILE"=DOWN10_quantile, "SCALE"=DOWN10_scale, "CYCLICLOESS"=DOWN10_cyclicloess))

venn(list("NONE"=UP100_none, "QUANTILE"=UP100_quantile, "SCALE"=UP100_scale, "CYCLICLOESS"=UP100_cyclicloess))
venn(list("NONE"=DOWN100_none, "QUANTILE"=DOWN100_quantile, "SCALE"=DOWN100_scale, "CYCLICLOESS"=DOWN100_cyclicloess))

venn(list("NONE"=UP1000_none, "QUANTILE"=UP1000_quantile, "SCALE"=UP1000_scale, "CYCLICLOESS"=UP1000_cyclicloess))
venn(list("NONE"=DOWN1000_none, "QUANTILE"=DOWN1000_quantile, "SCALE"=DOWN1000_scale, "CYCLICLOESS"=DOWN1000_cyclicloess))


venn(list("NONE"=UP_none, "QUANTILE"=UP_quantile))
venn(list("NONE"=UP_none, "SCALE"=UP_scale))
venn(list("NONE"=UP_none, "CYCLICLOESS"=UP_cyclicloess))
venn(list("QUANTILE"=UP_quantile, "SCALE"=UP_scale))
venn(list("QUANTILE"=UP_quantile, "CYCLICLOESS"=UP_cyclicloess))
venn(list("SCALE"=UP_scale, "CYCLICLOESS"=UP_cyclicloess))

venn(list("NONE"=DOWN_none, "QUANTILE"=DOWN_quantile))
venn(list("NONE"=DOWN_none, "SCALE"=DOWN_scale))
venn(list("NONE"=DOWN_none, "CYCLICLOESS"=DOWN_cyclicloess))
venn(list("QUANTILE"=DOWN_quantile, "SCALE"=DOWN_scale))
venn(list("QUANTILE"=DOWN_quantile, "CYCLICLOESS"=DOWN_cyclicloess))
venn(list("SCALE"=DOWN_scale, "CYCLICLOESS"=DOWN_cyclicloess))


#Figura PCA
PCA_none2 <- ggplot(data=pca.data_none2, aes(x=pca_none$x[,1], y=pca_none$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep="")) + theme_bw() + ggtitle("NONE") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")
PCA_quantile2 <- ggplot(data=pca.data2_quantile, aes(x=pca_quantile$x[,1], y=pca_quantile$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_quantile[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_quantile[2], "%", sep="")) + theme_bw() + ggtitle("QUANTILE") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")
PCA_scale2 <- ggplot(data=pca.data2_scale, aes(x=pca_scale$x[,1], y=pca_scale$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_scale[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_scale[2], "%", sep="")) + theme_bw() + ggtitle("SCALE") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")
PCA_cyclicloess2 <- ggplot(data=pca.data2_cyclicloess, aes(x=pca_cyclicloess$x[,1], y=pca_cyclicloess$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_cyclicloess[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_cyclicloess[2], "%", sep="")) + theme_bw() + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")

ext_leyenda_PCA <- ggplot(data=pca.data_none2, aes(x=pca_none$x[,1], y=pca_none$x[,2], label=Sample, color=Condition)) + geom_point() + xlab(paste("PC1 - ", pca.var.per_none[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per_none[2], "%", sep="")) + theme_bw() + ggtitle("NONE") + theme(legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12))
ext_leyenda_PCA <- ext_leyenda_PCA
leyenda_PCA <- get_legend(ext_leyenda_PCA)

#scale_fill_discrete(labels = c("Control", "Shock no séptico", "Shock séptico"))

PCA_grupo <- plot_grid(PCA_none2, NULL, PCA_quantile2, PCA_scale2, NULL, PCA_cyclicloess2, ncol = 3, nrow = 2, labels = c("A", "", "B", "C", "", "D"), rel_widths = c(1, 0.1, 1))
PCA_figura <- plot_grid(leyenda_PCA, PCA_grupo, ncol=1, rel_heights = c(1, 6))
PCA_figura


#Figura mapa de calor
heatmap_none2 <- pheatmap(df_none_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE, show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "None", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)
heatmap_quantile2 <- pheatmap(df_quantile_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE,show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Quantile", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)
heatmap_scale2 <- pheatmap(df_scale_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE,show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Scale", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)
heatmap_cyclicloess2 <- pheatmap(df_cyclicloess_short, cluster_rows = FALSE, cluster_cols = TRUE, border_color = FALSE,show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 8, annotation_col = annotation_df, main = "Cyclicloess", colorRampPalette(c("deepskyblue", "black", "yellow"))(25), scale = "row", annotation_legend = FALSE)

heatmap_grupo <- grid.arrange(heatmap_none2, NULL, heatmap_quantile2, heatmap_scale2, NULL, heatmap_cyclicloess2, ncol = 3, nrow = 2, labels = c("A", "", "B", "C", "", "D"), rel_widths = c(1, 0.1, 1))
heatmap_grupo


#Figura volcano
vp_none2 <- ggplot(data=toptablelimma_none, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("NONE") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")
vp_quantile2 <- ggplot(data=toptablelimma_quantile, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("QUANTILE") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")
vp_scale2 <- ggplot(data=toptablelimma_scale, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("SCALE") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")
vp_cyclicloess2 <- ggplot(data=toptablelimma_cyclicloess, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(legend.position = "none")

ext_leyenda_vp <- ggplot(data=toptablelimma_cyclicloess, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() + scale_color_manual(values=c("red", "black", "blue")) + geom_vline(xintercept=c(-0.5, 0.5), col="black") + geom_hline(yintercept=-log10(0.05), col="black") + ggtitle("CYCLICLOESS") + theme(plot.title = element_text(hjust=0.5)) + theme(legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12))
leyenda_vp <- get_legend(ext_leyenda_vp)

volcano_grupo <- plot_grid(vp_none2, NULL, vp_quantile2, vp_scale2, NULL, vp_cyclicloess2, ncol = 3, nrow = 2, labels = c("A", "", "B", "C", "", "D"), rel_widths = c(1, 0.1, 1))
volcanos_figura <- plot_grid(leyenda_vp, volcano_grupo, ncol=1, rel_heights = c(1, 6))
volcanos_figura









toptablelimma10_none <- toptablelimma_none[1:10,]
toptablelimma100_none <- toptablelimma_none[1:100,]
toptablelimma1000_none <- toptablelimma_none[1:1000,]
UP10_none <- row.names(toptablelimma10_none)[which(toptablelimma10_none$diffexpressed == "UP")]
UP100_none <- row.names(toptablelimma100_none)[which(toptablelimma100_none$diffexpressed == "UP")]
UP1000_none <- row.names(toptablelimma1000_none)[which(toptablelimma1000_none$diffexpressed == "UP")]
DOWN10_none <- row.names(toptablelimma10_none)[which(toptablelimma10_none$diffexpressed == "DOWN")]
DOWN100_none <- row.names(toptablelimma100_none)[which(toptablelimma100_none$diffexpressed == "DOWN")]
DOWN1000_none <- row.names(toptablelimma1000_none)[which(toptablelimma1000_none$diffexpressed == "DOWN")]

toptablelimma10_quantile <- toptablelimma_quantile[1:10,]
toptablelimma100_quantile <- toptablelimma_quantile[1:100,]
toptablelimma1000_quantile <- toptablelimma_quantile[1:1000,]
UP10_quantile <- row.names(toptablelimma10_quantile)[which(toptablelimma10_quantile$diffexpressed == "UP")]
UP100_quantile <- row.names(toptablelimma100_quantile)[which(toptablelimma100_quantile$diffexpressed == "UP")]
UP1000_quantile <- row.names(toptablelimma1000_quantile)[which(toptablelimma1000_quantile$diffexpressed == "UP")]
DOWN10_quantile <- row.names(toptablelimma10_quantile)[which(toptablelimma10_quantile$diffexpressed == "DOWN")]
DOWN100_quantile <- row.names(toptablelimma100_quantile)[which(toptablelimma100_quantile$diffexpressed == "DOWN")]
DOWN1000_quantile <- row.names(toptablelimma1000_quantile)[which(toptablelimma1000_quantile$diffexpressed == "DOWN")]

toptablelimma10_scale <- toptablelimma_scale[1:10,]
toptablelimma100_scale <- toptablelimma_scale[1:100,]
toptablelimma1000_scale <- toptablelimma_scale[1:1000,]
UP10_scale <- row.names(toptablelimma10_scale)[which(toptablelimma10_scale$diffexpressed == "UP")]
UP100_scale <- row.names(toptablelimma100_scale)[which(toptablelimma100_scale$diffexpressed == "UP")]
UP1000_scale <- row.names(toptablelimma1000_scale)[which(toptablelimma1000_scale$diffexpressed == "UP")]
DOWN10_scale <- row.names(toptablelimma10_scale)[which(toptablelimma10_scale$diffexpressed == "DOWN")]
DOWN100_scale <- row.names(toptablelimma100_scale)[which(toptablelimma100_scale$diffexpressed == "DOWN")]
DOWN1000_scale <- row.names(toptablelimma1000_scale)[which(toptablelimma1000_scale$diffexpressed == "DOWN")]

toptablelimma10_cyclicloess <- toptablelimma_cyclicloess[1:10,]
toptablelimma100_cyclicloess <- toptablelimma_cyclicloess[1:100,]
toptablelimma1000_cyclicloess <- toptablelimma_cyclicloess[1:1000,]
UP10_cyclicloess <- row.names(toptablelimma10_cyclicloess)[which(toptablelimma10_cyclicloess$diffexpressed == "UP")]
UP100_cyclicloess <- row.names(toptablelimma100_cyclicloess)[which(toptablelimma100_cyclicloess$diffexpressed == "UP")]
UP1000_cyclicloess <- row.names(toptablelimma1000_cyclicloess)[which(toptablelimma1000_cyclicloess$diffexpressed == "UP")]
DOWN10_cyclicloess <- row.names(toptablelimma10_cyclicloess)[which(toptablelimma10_cyclicloess$diffexpressed == "DOWN")]
DOWN100_cyclicloess <- row.names(toptablelimma100_cyclicloess)[which(toptablelimma100_cyclicloess$diffexpressed == "DOWN")]
DOWN1000_cyclicloess <- row.names(toptablelimma1000_cyclicloess)[which(toptablelimma1000_cyclicloess$diffexpressed == "DOWN")]






saveRDS(x, "x.rds")


venn(list("NONE"=UP_none, "QUANTILE"=UP_quantile))


ven <- venndetail(list("NONE"=UP_none, "QUANTILE"=UP_quantile, "SCALE"=UP_scale))
plot(ven)
plot(ven, type = "vennpie")



ven2 <- venn.diagram(list("NONE" = UP_none, "QUANTILE" = UP_quantile, "SCALE" = UP_scale), filename = "pruebaVenn.png", imagetype = "png", main="Diagrama de Venn", fill=c("lightgreen", "lightblue", "lightsalmon"), col=c("lightgreen", "lightblue", "lightsalmon"), cat.col=c("green", "blue", "salmon"))
  

venn.plot <- venn.diagram(list(UP_none, UP_quantile, UP_scale), NULL, imagetype="png", height = 480 , 
                          width = 480 , 
                          resolution = 600,
                          compression = "lzw",
                          lwd = 1,fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)), alpha=c(0.5,0.5,0.5), cex = 1,fontfamily = "sans", cat.col = c("#440154ff", '#21908dff', '#fde725ff'),cat.cex = 0.6,  cat.dist = c(0.055, 0.055, 0.085),cat.default.pos = "outer",cat.pos = c(-27, 27, 135),category.names=c("NONE","QUANTILE","SCALE"))



venn.plot <- venn.diagram(list(UP_none, UP_quantile, UP_scale), NULL, lwd = 1, fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)), alpha=c(0.5,0.5,0.5), cex = 1,fontfamily = "sans", cat.col = c("#440154ff", '#21908dff', '#fde725ff'),cat.cex = 0.6,  cat.dist = c(0.055, 0.055, 0.085),cat.default.pos = "outer",cat.pos = c(-27, 27, 135),category.names=c("NONE","QUANTILE","SCALE"))


venn.plot <- venn.diagram(list(UP_none, UP_quantile, UP_scale), NULL, lwd = 1, fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)), alpha=c(0.5,0.5,0.5), cex = 1, cat.col = c("#440154ff", '#21908dff', '#fde725ff'), cat.cex = 0.6,  cat.dist = c(0.055, 0.055, 0.085),cat.default.pos = "outer",cat.pos = c(-27, 27, 135),category.names=c("NONE","QUANTILE","SCALE"))


venn.plot_UP <- venn.diagram(list(UP_none, UP_quantile, UP_scale), NULL, lwd = 3, fill=c("lightgreen", "blue", "lightsalmon"), cat.default.pos = "text", category.names=c("NONE","QUANTILE","SCALE"))

venn.plot_DOWN <- venn.diagram(list(DOWN_none, DOWN_quantile, DOWN_scale), NULL, lwd = 3, fill=c("lightgreen", "blue", "lightsalmon"), cat.default.pos = "text", category.names=c("NONE","QUANTILE","SCALE"))

grid.newpage()
grid.draw(venn.plot_UP)

grid.newpage()
grid.draw(venn.plot_DOWN)

grid.newpage()
draw.triple.venn(list("NONE" = UP_none, "QUANTILE" = UP_quantile, "SCALE" = UP_scale), n12 = length(intersect(UP_none, UP_quantile)), n23 = length(intersect(UP_quantile, UP_scale)), n13 = length(intersect(UP_none, UP_scale)), n123 = length(intersect(intersect(UP_none, UP_quantile), UP_scale)))


  vlist,     
             filename="Venn_3way_more.png",
             imagetype="png",
             main="Venn diagram",
             sub="3-way",
             main.col="red",
             fill=c("lightgreen", "lightblue", "lightsalmon"),
             col=c("lightgreen", "lightblue", "lightsalmon"),
             cat.col=c("green", "blue", "salmon"))



intersect(UP_none, UP_quantile)
intersect(UP_none, UP_scale)
intersect(UP_quantile, UP_scale)
intersect(intersect(UP_none, UP_quantile), UP_scale)

intersect(DOWN_none, DOWN_quantile)
intersect(DOWN_none, DOWN_scale)
intersect(DOWN_quantile, DOWN_scale)
intersect(intersect(DOWN_none, DOWN_quantile), DOWN_scale)




names_none <- rownames(toptablelimma_none)
toptablelimma_noneVD <- cbind(toptablelimma_none, names_none)
colnames(toptablelimma_noneVD)[7] <- "diffexpressed_none"
toptablelimma_noneVD <- toptablelimma_noneVD[, c(1,7)]

names_quantile <- rownames(toptablelimma_quantile)
toptablelimma_quantileVD <- cbind(toptablelimma_quantile, names_quantile)
colnames(toptablelimma_quantileVD)[7] <- "diffexpressed_quantile"

toptablelimma_noneVD <- toptablelimma_none[, c(1,7)]
toptablelimma_quantileVD <- toptablelimma_quantile[, c(1,7)]

pruebaVD <- cbind(toptablelimma_noneVD, toptablelimma_quantileVD)

toptablelimma_noneVD[, c(1,7)]

VD <- list(toptablelimma_none$diffexpressed == "UP", toptablelimma_quantile$diffexpressed == "UP")
VD

none_UP 




length(UP_quantile)

length(toptablelimma_none$diffexpressed == "DOWN")



#Diagrama de Venn
dge2 <- CvsSS_dge[CvsSS_dge$table$FDR <= 0.05 & abs(CvsSS_dge$table$logFC) >= 1,]
dim(dge2)
dge_up <- rownames(dge2[dge2$table$logFC>=1,])
dge_down <- rownames(dge2[dge2$table$logFC<=-1,])
length(dge_up)
length(dge_down)
venn(list("Upregulated"=dge_up, "Downregulated"=dge_down))

A_class <- decideTestsDGE(CvsSS)
summary(A_class)

plotMD(CvsSS, status=A_class, main="SS vs NSS", cex=0.2)







toptablelimma_none2 <- rownames_to_column(var = "Probe")



venn.plot <- venn.diagram(list(UP_none, UP_quantile, UP_scale), NULL,
                          height = 480 , 
                          width = 480 , 
                          resolution = 600,
                          compression = "lzw",
                          lwd = 1,fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)), alpha=c(0.5,0.5,0.5), cex = 1,fontfamily = "sans", cat.col = c("#440154ff", '#21908dff', '#fde725ff'), cat.cex = 0.6, cat.dist = c(0.055, 0.055, 0.085), cat.default.pos = "outer",cat.pos = c(-27, 27, 135),category.names=c("NONE","QUANTILE","SCALE"))


tiff("./Desktop/Venn_Diagram_1.tiff")
grid.draw(venn.plot)
dev.off()

ven <- venndetail(list(Cortex = T2DM$Cortex$Entrez, SCN = T2DM$SCN$Entrez, Glom = T2DM$Glom$Entrez))




venn.plot <- draw.pairwise.venn(toptablelimma_none$diffexpressed == "UP", toptablelimma_quantile$diffexpressed == "UP", cross.area = TRUE)


  area1 = 100,
  area2 = 70,
  cross.area = 68,
  category = c("First", "Second"),
  fill = c("blue", "red"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(285, 105),
  cat.dist = 0.09,
  cat.just = list(c(-1, -1), c(1, 1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed"
);
grid.draw(venn.plot);
grid.newpage()



draw.pairwise.venn()

UP <- cbind(UP_none, UP_quantile)

vd <- venn.diagram(list(UP_none, UP_quantile), fill = c("red", "green"), alpha = c(0.5, 0.5), filename=NULL)

vd <- venn.diagram(list(UP_none, UP_quantile), fill = c("red", "green"), alpha = c(0.5, 0.5), filename=NULL)

vd <- venn.diagram(list(UP=toptablelimma_none$diffexpressed, UP=toptablelimma_quantile$diffexpressed), fill = c("red", "green"), alpha = c(0.5, 0.5), filename=NULL)
vd




# limma-Voom
dgeLV <- DGEList(data_final)
normLV <- calcNormFactors(dgeLV)
plotMDS(normLV, col = as.numeric(diagnostico))
plotMDS(normLV)
V <- voom(normLV, diseño, plot=TRUE)
fitlvoom <- lmFit(V, diseño)
fitlvoom <- contrasts.fit(fitlvoom, contraste)
fitlvoom <- eBayes(fitlvoom)
toptableLV <- topTable(fitlvoom, sort.by = "P", n = Inf)
toptableLV
resultsLV <- decideTests(fitlvoom)
vennDiagram(resultsLV)
length(which(toptableLV$adj.P.Val < 0.05))
siglV <- subset(toptableLV,adj.P.Val<0.05)
siglV_up <- rownames(subset(sigLV,logFC>0))
siglV_dn <- rownames(subset(sigVL,logFC<0))
length(siglV_up)
length(siglV_dn)



# creating venn diagram for four sets
ggvenn(list(siglimma_up, sigVL_up), c("siglimma_up", "sigVL_up"), show_percentage=TRUE, fill_color=c("red","orange"))


grid.newpage()
venn.diagram(list(siglimma_up, sigVL_up), filename = NULL, category.names = c("Set 1" , "Set 2"), fill = c("#E69F00", "#009E73"))
draw.single.venn(area = length(siglimma_up))



coloresvenn <- c("#6b7fff","#c3db0f")
venn.diagram(list(siglimma_up, sigVL_up), filename = tempfile(pattern = "pruebavenn", fileext = ".png"), category.names=c("Limma_up", "Voom/Limma_up"), col="black", fill= coloresvenn, cat.col=coloresvenn)


venn(list("Limma_up"=siglimma_up, "Voom/Limma_up"=sigVL_up))
venn(list("Limma_dn"=siglimma_dn, "Voom/Limma_dn"=sigVL_dn))





# EdgeR LRT
data_EdgeR <- data_final[ ,-c(115:129)]
metadatos_EdgeR <- metadatos[-c(115:129),]
grupo_EdgeR <- paste(metadatos_EdgeR$diagnosis, sep=";")
grupo_EdgeR <- factor(grupo_EdgeR)
table(grupo_EdgeR)
diseño_EdgeR <- model.matrix(~0+grupo_EdgeR)
diseño_EdgeR
colnames(diseño_EdgeR) <- levels(grupo_EdgeR)
rownames(diseño_EdgeR) <- rownames(metadatos_EdgeR)
diagnostico_EdgeR <- paste(metadatos_EdgeR$diagnosis)
diagnostico_EdgeR


dgeElrt <- DGEList(data_EdgeR)
normElrt <- calcNormFactors(dgeElrt)
dispElrt <- estimateDisp(normElrt, diseño_EdgeR, robust=TRUE, prior.df=1)
fitElrt <- glmFit(dispElrt, diseño_EdgeR)
fitElrt <- glmLRT(fitElrt)

dge<-as.data.frame(topTags(fitElrt,n=Inf))
dge$dispersion<-fitElrt$dispersion
dge<-merge(dge,fitElrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-dge[order(dge$PValue),]
head(dge,10)
dge_edger <- dge
sig <- subset(dge_edger,FDR<0.05)
dge_edger_up <- rownames(subset(dge_edger,logFC>0))
dge_edger_dn <- rownames(subset(dge_edger,logFC<0))
length(dge_edger_up)
length(dge_edger_dn)


# EdgeR QL
y <- DGEList(counts=data_final)
y <- calcNormFactors(y)
y <- estimateDisp(y, design=diseño, robust=TRUE,prior.df=1)
fitql <- glmQLFit(y, design=diseño)
lrtql <- glmQLFTest(fitql)
dgeql <- as.data.frame(topTags(lrtql,n=Inf))
dgeql$dispersion <- lrtql$dispersion
dgeql <- merge(dgeql,lrtql$fitted.values,by='row.names')
rownames(dgeql)=dgeql$Row.names
dgeql$Row.names=NULL
dgeql <- dgeql[order(dgeql$PValue),]
head(dgeql,10)
dge_edgerql <- dgeql
sigql <- subset(dge_edgerql, FDR<0.05)
dge_edgerql_up <- rownames(subset(sigql,logFC>0))
dge_edgerql_dn <- rownames(subset(sigql,logFC<0))
length(dge_edgerql_up)
length(dge_edgerql_dn)









volcanoplot(fitlimma, names = fitlimma$genes$Gene.symbol, highlight = 5)




#DESeq2
data2 <- makeSummarizedExperimentFromDataFrame(data_final, diseño, ignore.strand = TRUE)




dds <- DESeqDataSet(data2, design = diseño)
res <- DESeq(dds)

z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
head(dge,10)
dge_deseq2 <- dge
sig <- subset(dge,padj<0.05)
dge_deseq2_up <- rownames(subset(sig,log2FoldChange>0))
dge_deseq2_dn <- rownames(subset(sig,log2FoldChange<0))
length(dge_deseq2_up)
length(dge_deseq2_dn)




#Creación de la lista de DEG
deg <- DGEList(data_final)
dim(deg)
class(deg)
names(deg)


#Adición de la información de cada una de las muestras
grupo 
deg$samples$group <- grupo
deg$samples


#Normalización de los datos mediante el método TMM que de manera predeterminada utiliza la función calcNormFactors
deg <- calcNormFactors(deg)
deg$samples


#Representación con los valores promedio de la expresión logarítmica (valores M) para cada una de las muestras. Proviene de RNAseq
cpms = cpm(deg,log=TRUE)
boxplot(cpms, xlab="muestras", ylab="Log2 CPM",las=2)
abline(h=median(cpms),col="red")
title("Tamaños librerías (logCPMs)")


#Representación de la relación media-varianza de los datos
deg <- estimateDisp(deg, diseño, robust=TRUE)
#No acaba de funcionar
plotMeanVar(deg, show.tagwise.vars=TRUE, NBline=TRUE)


#Ajuste binomio negativo GLM para cada gen
fit <- glmQLFit(deg, diseño, robust=TRUE)
dim(fit)

#topTable de la comparación entre control y shock séptico
contraste <- makeContrasts(control-SS, levels = grupo)
fit2 <- lmFit(data_final, diseño)
fit3 <- contrasts.fit(fit2, contraste)
fitted.eBayes <- eBayes(fit3)
tab <- topTable(fitted.eBayes, n=Inf, coef=1)




#Representación de la dispersión de cada gen frente a la expresion media
plotQLDisp(fit, ylab = "dispersión QL", ylim=c(0,4), xlim=c(-3,12))


#Prueba de la hipótesis nula entre control y shock séptico
CvsSS <- glmQLFTest(fit, contrast=c(-1,0,1))
head(CvsSS)


#Prueba de la hipótesis nula entre control y shock no séptico
CvsNSS <- glmQLFTest(fit, contrast=c(-1,1,0))
head(CvsNSS)


#Prueba de la hipótesis nula entre shock no séptico y shock séptico
NSSvsSS <- glmQLFTest(fit, contrast=c(0,1,-1))
head(NSSvsSS)


#A partir de este punto se trabaja sólo con la comparación entre control y shock séptico para simplificar

#Obtención de los genes diferencialmente expresados en la comparación establecida. Se puede elegir el valor del p-valor de corte
CvsSS_deg <- topTags(CvsSS, n = Inf, sort.by = "PValue")
head(CvsSS_deg)


#Cambios en los histogramas del p-valor y de FDR. Comprobar si la gráfica de FDR es correcta
hist(CvsSS_deg$table$PValue, breaks=50, xlab ="p-value", main = "Histogram of p-values")
hist(CvsSS_deg$table$FDR, breaks=50, xlab ="p-value", main = "Histogram of FDR")



#Diagrama de Venn. Valor logFC modificado a muy bajo para que funcione
dge2 <- CvsSS_deg[CvsSS_deg$table$FDR <= 0.05 & abs(CvsSS_deg$table$logFC) >= 0.5,]
dim(dge2)
dge_up <- rownames(dge2[dge2$table$logFC>=0.5,])
dge_down <- rownames(dge2[dge2$table$logFC<=-0.5,])
length(dge_up)
length(dge_down)
venn(list("Upregulated"=dge_up, "Downregulated"=dge_down))

A_class <- decideTestsDGE(CvsSS)
summary(A_class)
plotMD(CvsSS, status=A_class, main="SS vs NSS", cex=0.2)


#VOLCANO PLOT (NO ES OPCIÓN GGPLOT). Algo sucede porque no la gran mayoría de los genes están en valor 0 de -log10FDR
vol <- as.data.frame(cbind(CvsSS_deg$table$logFC, -log10(CvsSS_deg$table$FDR)))
colnames(vol) <- c("logFC", "-log10FDR")
rownames(vol) <- rownames(CvsSS_deg)
#Selección de etiquetas
n <- head(vol[order(-vol$`-log10FDR`),], 10)
#Rojo para los DEG
colores <- ifelse(CvsSS_deg$table$FDR <=0.05 & CvsSS_deg$table$logFC >=1 | CvsSS_deg$table$FDR <= 0.05 & CvsSS_deg$table$logFC <=-1 , "red", "black")
#Plot
volcan <- plot(vol, pch = 19, col = colores, cex = 0.2)
text(n, labels=rownames(n), cex=0.5)


#Prueba de heatmap con las cien primeras filas (genes) para evitar perder mucho tiempo en cálculo
cienfilas <- as.matrix(data_final_none[1:100,])
colores <- colorRampPalette(brewer.pal(11, "RdBu"))
#colores <- colorRampPalette(c("red", "green"))
#coloresanotacion <- list(diseño = c(SS = "chartreuse4", NSS = "burlywood3"))
#Añadir a heatmap annotation_colors = coloresanotacion, pero no acaba de funcionar bien.
heatmap.2(cienfilas, col = colores, main = "Heatmap", trace = "none", margins = c(10,12), cexRow = 0.7, Rowv = TRUE, Colv = TRUE,  key = FALSE, xlab = "gen", ylab = "paciente")



prueba_none <- as.matrix(data_final_none)
colores <- colorRampPalette(brewer.pal(11, "RdBu"))
heatmap(prueba_none)
heatmap.2(prueba_none, col = colores, main = "NONE", trace = "none", margins = c(10,12), cexRow = 0.7, Rowv = FALSE, Colv = FALSE,  key = FALSE)




#Establecimiento de la comparación
contraste <- makeContrasts(control-SS, levels = grupo)

#como alternativa se puede probar BAYES y ver si hay diferencias (minuto 24 video de ayuda)
CvsSS <- glmQLFTest(fit, contrast = contraste)
dim(CvsSS)

#Alternativa a glmQLFTest, no hay grandes cambios
#CvsSS2 <- glmLRT(fit, contrast = contraste)
#CvsSS_dge2 <- topTags(CvsSS2, n=Inf)


#TopTable 
CvsSS_dge <- topTags(CvsSS, n=Inf)
dim(CvsSS)

#Representación de la toptag (y en orden ascendente y descendente)
head(CvsSS_dge)
head(CvsSS_dge$table[order(-CvsSS_dge$table$F),])
head(CvsSS_dge$table[order(CvsSS_dge$table$F),])
head(CvsSS_dge$table[order(CvsSS_dge$table$logFC),])
head(CvsSS_dge$table[CvsSS_dge$table$logFC==0,])

#Diagrama de Venn
dge2 <- CvsSS_dge[CvsSS_dge$table$FDR <= 0.05 & abs(CvsSS_dge$table$logFC) >= 1,]
dim(dge2)
dge_up <- rownames(dge2[dge2$table$logFC>=1,])
dge_down <- rownames(dge2[dge2$table$logFC<=-1,])
length(dge_up)
length(dge_down)
venn(list("Upregulated"=dge_up, "Downregulated"=dge_down))

A_class <- decideTestsDGE(CvsSS)
summary(A_class)

plotMD(CvsSS, status=A_class, main="SS vs NSS", cex=0.2)


#VOLCANO PLOT (NO ES OPCIÓN GGPLOT)
vol <- as.data.frame(cbind(CvsSS_dge$table$logFC, -log10(CvsSS_dge$table$FDR)))
colnames(vol) <- c("logFC", "-log10FDR")
rownames(vol) <- rownames(CvsSS_dge)
#Selección de etiquetas
n <- head(vol[order(-vol$`-log10FDR`),], 10)
#Rojo para los DEG
colores <- ifelse(CvsSS_dge$table$FDR <=0.05 & CvsSS_dge$table$logFC >=1 | CvsSS_dge$table$FDR <= 0.05 & CvsSS_dge$table$logFC <=-1 , "red", "black")
#Plot
volcan <- plot(vol, pch = 19, col = colores, cex = 0.2)
text(n, labels=rownames(n), cex=0.5)



#Prueba de heatmap con las cien primeras filas (genes) para evitar perder mucho tiempo en cálculo
ciengenes <- as.matrix(data_final[1:100,])
colores <- colorRampPalette(brewer.pal(11, "RdBu"))
#colores <- colorRampPalette(c("red", "green"))
#coloresanotacion <- list(diseño = c(SS = "chartreuse4", NSS = "burlywood3"))
#Añadir a heatmap annotation_colors = coloresanotacion, pero no acaba de funcionar bien.
heatmap.2(ciengenes, col = colores, main = "Prueba", trace = "none", margins = c(10,12), cexRow = 0.7, Rowv = TRUE, Colv = TRUE,  key = TRUE, xlab = "gen", ylab = "paciente")







CvsSS_deg$table$exp <- ifelse(CvsSS_deg$table$logFC > 0.5 & CvsSS_deg$table$PValue <= 0.05, "up_volcano", ifelse(CvsSS_deg$table$logFC < -0.5 & CvsSS_deg$table$PValue <= 0.05, "down_volcano", ifelse(CvsSS_deg$table$logFC [-0.5:0.5] & CvsSS_deg$table$PValue > 0.05, "no_volcano")))
pruebacolores <- case_when(CvsSS_deg$table$logFC >= 0.5 & CvsSS_deg$table$PValue <= 0.05 ~ "up_volcano", CvsSS_deg$table$logFC <= 0.5 & CvsSS_deg$table$PValue <= 0.05 ~ "down_volcano", TRUE ~ "ns_volcano")
colores2 <- c("red", "green", "black")
names(colores2) <- c("down_volcano", "up_volcano", "no_volcano")
volcano <- ggplot(data=CvsSS$table, aes(x=logFC, y=-log10(PValue), colour = CvsSS_deg$table$exp)) + geom_point(alpha=0.4, size=1.75) + geom_vline(xintercept=c(-0.5, 0.5), col="grey", linetype = "dashed") + geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed")+ theme_minimal()
volcano












#Prueba de heatmap con las cien primeras filas (genes) para evitar perder mucho tiempo en cálculo
cienfilas <- as.matrix(data_final_none[1:100,])
colores <- colorRampPalette(brewer.pal(11, "RdBu"))
#colores <- colorRampPalette(c("red", "green"))
#coloresanotacion <- list(diseño = c(SS = "chartreuse4", NSS = "burlywood3"))
#Añadir a heatmap annotation_colors = coloresanotacion, pero no acaba de funcionar bien.
heatmap.2(cienfilas, col = colores, main = "Heatmap", trace = "none", margins = c(10,12), cexRow = 0.7, Rowv = TRUE, Colv = TRUE,  key = FALSE, xlab = "gen", ylab = "paciente")





heatmap(prueba_none)


heatmap.2(prueba_none, col = colores, main = "NONE", trace = "none", margins = c(10,12), cexRow = 0.7, Rowv = TRUE, Colv = TRUE,  key = TRUE, labCol = FALSE)
heatmap.2(prueba_none, col = colores, main = "NONE", dendrogram = "col", trace = "none", margins = c(10,12), cexRow = 0.7, key = TRUE, key.title = "KEY", labCol = FALSE) 


heatmap.2(toptablelimma_none, col = colores, main = "NONE", trace = "none", margins = c(10,12), cexRow = 0.7, Rowv = TRUE, Colv = TRUE,  key = FALSE)


pheatmap(prueba_none, clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan", show_colnames=FALSE)


hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)
pheatmap(dists, col = rev(hmcol), clustering_distance_rows = "manhattan",show_colnames=F,fontsize_row=3,fontsize=3,clustering_distance_cols = "manhattan", treeheight_col = 0)
f<-factor(targets$Condition, levels=unique(targets$Condition))
design<-model.matrix(~0+f)
colnames(design)<-levels(f)
















#Anotación de base de datos
anotados <- data_agreg3 %>% rownames_to_column(var = "ProbeName") %>% inner_join(., featuredata, by = "ID")



data2 <- column_to_rownames(data, var = "ProbeName")

#No da error, hay que comprobar los nombres de las variables MIRAR EN SCRIPT DESARROLLADO




anotadosagregados <- aggregate(e.raw_norm, sum)


anotados <- data %>% rownames_to_column(var = "ProbeName") %>% inner_join(., featuredata, by = "ID")

data_agreg <- aggregate(data, list (by = featuredata, FUN = mean)
data_agreg_col <- data_agreg %>% remove_rownames() %>% column_to_rownames(var = "Group.1")
data_agreg_col <- data_agreg_col[,-c(1)]



#Control de calidad NO FUNCIONA
#agQuality(fnames = targets, organism = c("Mm", "Hs"), compBoxplot = TRUE, reference = NULL, controlMatrix = agcontrolCode, controlId = c("ProbeName"), output = FALSE, resdir = ".", dev= "png", DEBUG = FALSE)


boxplot(featuredata)         
