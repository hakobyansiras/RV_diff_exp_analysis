library(DESeq2)
library(data.table)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(psf)
library(reshape2)
library(ggpubr)
library(googlesheets4)
library(fgsea)
library(UpSetR)

### Reading count matrix and sample information
rv_counts <- read.delim(file = "rv_raw_counts_3prime_500.tsv", sep = "\t")
coldata <- read.delim(file = "rv_coldata.tsv", sep = "\t", stringsAsFactors = F)

## Low count filtering
rv_counts <- rv_counts[,coldata$sample_id]
rv_counts_filtered <- rv_counts[which(rowSums(rv_counts >= 20) >= 56),]


### Assigning levels to experimental groups
coldata$gender <- factor(coldata$gender, levels = c("male", "female"))
coldata$time_collect <- factor(coldata$time_collect, levels = c("young", "old"))
coldata$Field <- factor(coldata$Field, levels = c("Sham", "Gamma", "Mix"))
coldata$group_new <- factor(coldata$group_new, levels = c("Sham male 440d", "Sham female 440d", 
                                                          "Sham male 660d", "Sham female 550d", 
                                                          "Gamma male 440d", "Gamma female 440d", 
                                                          "simGCRsim male 440d",  "simGCRsim female 440d",
                                                          "Gamma male 660d", "Gamma female 550d", 
                                                          "simGCRsim male 660d",  "simGCRsim female 550d"
))

### Plotting RIN values for each experimental group
layout(matrix(1:2, ncol = 2))
par(mar = c(10.5, 4, 4, 2))
boxplot(coldata$RIN ~ coldata$group_new, las = 2, xlab = "", ylab = "RIN", main = "RV rad RIN")

barplot(
  sapply(levels(coldata$group_new), function(x) {
    rins <- coldata[which(coldata$group_new == x),"RIN"]
    if(length(rins) < 5) {
      c(rins, rep(NA, 5-length(rins)))
    } else {
      rins
    }
  }), beside = T, col = "grey", las = 2, main = "RV rad RIN"
)

### Merging age groups togather
coldata$group_age_combined <- paste(stringr::str_to_title(coldata$gender), gsub("Mix", "simGCRsim", coldata$Field))
coldata$group_age_combined <- factor(coldata$group_age_combined, levels = c("Female Sham", "Female Gamma", "Female simGCRsim", "Male Sham", "Male Gamma", "Male simGCRsim"))
coldata <- coldata[order(coldata$group_age_combined),]


### Running DESeq with comined age groups
dds <- DESeqDataSetFromMatrix(countData = rv_counts_filtered[,coldata$sample_id], 
                              colData = coldata, 
                              design = ~RIN+group)

dds <- DESeq(dds)


### PCA plot on vst normalized expresssion data
palette_colors <- RColorBrewer::brewer.pal(6, "Set1")

plotPCA(vst(dds, blind=FALSE), intgroup="group_age_combined") + 
  scale_color_manual(values=palette_colors) +  # Apply color palette
  theme_minimal(base_size = 14) +  # Improve aesthetics
  labs(title="PCA Plot") +
  theme(
    legend.position="right",
    legend.title = element_blank(),
    panel.grid.major = element_line(color="grey90"),
    panel.grid.minor = element_blank()
  )

ggsave("Figure1.png",
       units="in", width=10/1.5, height=7/1.5, dpi=320, device = "png", bg = "white")

### Pefroming ANOVA test on PCs
right_ventrical_rad_norm_counts <- assay(vst(dds, blind=F))

pca_result <- prcomp(t(right_ventrical_rad_norm_counts), scale. = TRUE)

rad_groups <- setNames(nm = coldata$sample_id, object = as.character(coldata$Field))[rownames(pca_result$x)]
age_groups <- setNames(nm = coldata$sample_id, object = as.character(coldata$time_collect))[rownames(pca_result$x)]
sex_groups <- setNames(nm = coldata$sample_id, object = as.character(coldata$gender))[rownames(pca_result$x)]
RIN <- setNames(nm = coldata$sample_id, object = as.character(coldata$RIN))[rownames(pca_result$x)]
QPCR_mol <- setNames(nm = coldata$sample_id, object = as.character(coldata$QPCR_mol))[rownames(pca_result$x)]
concentration <- setNames(nm = coldata$sample_id, object = as.character(coldata$concentration))[rownames(pca_result$x)]

aov_res <- lapply(c("rad_groups", "age_groups", "sex_groups", "RIN", "QPCR_mol", "concentration"), function(z) {
  
  sort(sapply(colnames(pca_result$x), function(y) {
    
    aov_res <- aov(pca_result$x[,y] ~ get(z))
    
    summary(aov_res)[[1]][["Pr(>F)"]][1]
  }))
  
})

names(aov_res) <- c("rad_groups", "age_groups", "sex_groups", "RIN", "QPCR_mol", "concentration")

print(aov_res)



### Performing differential expression analysis for each gender separately

## Females
rv_rad_coldata_females <- coldata[which(as.character(coldata$gender) == "female"),]

rv_rad_coldata_females$group <- factor(rv_rad_coldata_females$group, levels = c("Sham_female_young",
                                                                                "Sham_female_old",
                                                                                "Mix_female_young",
                                                                                "Mix_female_old",
                                                                                "Gamma_female_young",
                                                                                "Gamma_female_old"
))


rv_females_dds <- DESeqDataSetFromMatrix(countData = rv_counts_filtered[,rv_rad_coldata_females$sample_id], 
                                         colData = rv_rad_coldata_females, 
                                         design = ~Field+time_collect+RIN+QPCR_mol+concentration)
rv_females_dds <- DESeq(rv_females_dds)

resultsNames(rv_females_dds)

female_RIN_degs <- as.data.frame(results(rv_females_dds, name = "RIN"))
female_RIN_degs <- female_RIN_degs[which(female_RIN_degs$padj < 0.05),]


### Males
rv_rad_coldata_males <- coldata[which(as.character(coldata$gender) == "male"),]
rv_rad_coldata_males$group <- factor(rv_rad_coldata_males$group, levels = c("Sham_male_young", "Sham_male_old", "Mix_male_young", "Mix_male_old", "Gamma_male_young", "Gamma_male_old"))




rv_males_dds <- DESeqDataSetFromMatrix(countData = rv_counts_filtered[,rv_rad_coldata_males$sample_id], 
                                       colData = rv_rad_coldata_males, 
                                       design = ~Field+time_collect+RIN+QPCR_mol+concentration)
rv_males_dds <- DESeq(rv_males_dds)

resultsNames(rv_males_dds)


load("mouse_gene_conversions.RData")


male_RIN_degs <- as.data.frame(results(rv_males_dds, name = "RIN"))
male_RIN_degs <- male_RIN_degs[which(male_RIN_degs$padj < 0.05),]


male_mix_degs <- as.data.frame(results(rv_males_dds, contrast=c("Field","Mix","Sham")))
male_mix_degs <- male_mix_degs[which(male_mix_degs$padj < 0.05),]
male_mix_degs <- male_mix_degs[order(male_mix_degs$padj),]

male_mix_degs$gene_name <- gene_id_to_symbol[rownames(male_mix_degs)]
male_mix_degs$gene_type <- gene_id_to_type[rownames(male_mix_degs)]
male_mix_degs$gene_description <- gene_id_to_description[rownames(male_mix_degs)]
male_mix_degs$chr <- gene_id_to_chromosome[rownames(male_mix_degs)]

male_mix_degs <- male_mix_degs[setdiff(rownames(male_mix_degs), rownames(male_RIN_degs)),]



### Over represenation analysis of Male simGCRsim IR group DEGs
library(enrichR)

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", 
         "KEGG_2019_Mouse")

male_mix_ORA <- enrichr(gene_id_to_symbol[setdiff(rownames(male_mix_degs), rownames(male_RIN_degs))], dbs) 

sapply(male_mix_ORA, function(x) {sum(x$Adjusted.P.value < 0.05)})


### Extracting DEGs from DESeq results
rv_rad_deg_tables <- c(
  lapply(list(Female_Gamma = as.data.frame(results(rv_females_dds, contrast=c("Field","Gamma","Sham"))),
              Female_Mix = as.data.frame(results(rv_females_dds, contrast=c("Field","Mix","Sham")))), function(x) {
                x$gene_name <- gene_id_to_symbol[rownames(x)]
                x$gene_type <- gene_id_to_type[rownames(x)]
                x$gene_description <- gene_id_to_description[rownames(x)]
                x$chr <- gene_id_to_chromosome[rownames(x)]
                x$RIN_overlap <- ifelse(rownames(x) %in% rownames(female_RIN_degs), "Yes", "")
                x[order(x$padj),]
              }),
  lapply(list(Male_Gamma = as.data.frame(results(rv_males_dds, contrast=c("Field","Gamma","Sham"))),
              Male_Mix = as.data.frame(results(rv_males_dds, contrast=c("Field","Mix","Sham")))), function(x) {
                x$gene_name <- gene_id_to_symbol[rownames(x)]
                x$gene_type <- gene_id_to_type[rownames(x)]
                x$gene_description <- gene_id_to_description[rownames(x)]
                x$chr <- gene_id_to_chromosome[rownames(x)]
                x$RIN_overlap <- ifelse(rownames(x) %in% rownames(male_RIN_degs), "Yes", "")
                x[order(x$padj),]
              })
)

### Plotting volcano plots for RV comparison groups
y_lim <- c(0, max(unlist(lapply(rv_rad_deg_tables, function(x) {-log10(x$padj)})), na.rm = T))
x_lim <- range(unlist(lapply(rv_rad_deg_tables, function(x) {x$log2FoldChange})), na.rm = T)

rv_volcanos <- lapply(names(rv_rad_deg_tables), function(x) {
  up_num <- sum(rv_rad_deg_tables[[x]]$log2FoldChange[which(rv_rad_deg_tables[[x]]$padj < 0.05)] > 0)
  down_num <- sum(rv_rad_deg_tables[[x]]$log2FoldChange[which(rv_rad_deg_tables[[x]]$padj < 0.05)] < 0)
  
  EnhancedVolcano(rv_rad_deg_tables[[x]],
                  lab = rv_rad_deg_tables[[x]]$gene_name,
                  ylim = y_lim, xlim = x_lim,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = x, # if(grep(x, names(diff_exps_deseq_full)[c(6,5,8,7)]) > 2) {NULL} else {titles[grep(x, names(diff_exps_deseq_full)[c(6,5,8,7)])]},
                  subtitle = paste0("Up ", up_num, " Down ", down_num),
                  subtitleLabSize = 16,
                  caption = NULL,
                  pCutoff = 0.05,
                  FCcutoff = 0.5,
                  drawConnectors = TRUE,
                  legendLabSize = 10,
                  captionLabSize = 8,
                  pointSize = 2.0,
                  labSize = 5.0) + theme(plot.title = element_text(hjust = 0.5))
})


patchwork::wrap_plots(rv_volcanos)




rv_rad_filtered_degs <- lapply(rv_rad_deg_tables, function(x) {
  p_val_filtered <- x[which(x$padj < 0.05),]
  p_val_filtered[which(p_val_filtered$RIN_overlap == ""),]
})

rv_rad_filtered_degs <- lapply(rv_rad_filtered_degs, function(x) {
  up <- x[which(x$log2FoldChange > 0),]
  down <- x[which(x$log2FoldChange < 0),]
  
  deg_table <- rbind(up[order(up$log2FoldChange, decreasing = T),],
                     down[order(down$log2FoldChange),])
  
  deg_table$log2FoldChange <- round(deg_table$log2FoldChange, digits = 3)
  
  deg_table$gene_description <- gsub("\\[.*?\\]", "", deg_table$gene_description)
  
  deg_table[,c("gene_type", "gene_name", "gene_description", "log2FoldChange", "padj")]
})






### Expression heatmap for simGCRsim IR male DEGs
rv_rad_selected_heatmap_genes <- unique(unlist(sapply(rv_rad_filtered_degs, function(x) {
  if(nrow(x) > 1) {
    fc_sorted_exp <- x[order(abs(x$log2FoldChange), decreasing = T),]
    rownames(fc_sorted_exp) # [1:min(nrow(fc_sorted_exp), 50)]
  }
})))

rv_rad_selected_heatmap_genes <- unique(unlist(sapply(rv_rad_filtered_degs, rownames)))

library(ComplexHeatmap)
ann <- HeatmapAnnotation(df = data.frame(Group = coldata$group_age_combined, row.names = coldata$sample_id), 
                         col = list(Group = setNames(object = palette_colors, nm = levels(coldata$group_age_combined))[as.character(coldata$group_age_combined)]))

png("Figure3.png",width=10/1.5,height=7/1.5,units="in",res=320, bg = "white")

Heatmap(assay(vst(dds, blind=FALSE))[rv_rad_selected_heatmap_genes,]/rowMeans(assay(vst(dds, blind=FALSE))[rv_rad_selected_heatmap_genes,]), top_annotation = ann,
        show_row_names = F,
        cluster_columns = F, name = "Norm counts",  heatmap_legend_param = list(labels_gp = gpar(fontsize = 12)) )

dev.off()




### Differential expression analysis of left ventircel IR groups to compare with right ventircle IR group DEGs

## Importing left ventircle raw counts and sample infromation
left_ventrical_raw_counts <- read.delim(file = "lv_raw_counts.tsv", sep = "\t")


lv_coldata <- read.delim(file = "lv_coldata.tsv", sep = "\t", stringsAsFactors = F)
lv_coldata$time_collect[which(lv_coldata$time_collect != "16 months")] <- "old"
lv_coldata$time_collect[which(lv_coldata$time_collect == "16 months")] <- "young"
lv_coldata$group <- paste(lv_coldata$Field, coldata$gender, coldata$time_collect, sep = "_")
lv_coldata <- lv_coldata[which(lv_coldata$Diet == "Normal"),]

lv_coldata$gender <- factor(lv_coldata$gender, levels = c("male", "female"))
lv_coldata$time_collect <- factor(lv_coldata$time_collect, levels = c("young", "old"))
lv_coldata$Field <- factor(lv_coldata$Field, levels = c("Sham", "Gamma", "Mix"))
lv_coldata$group <- factor(lv_coldata$group, levels = c("Sham_male_young", "Sham_female_young", 
                                                        "Sham_male_old", "Sham_female_old", 
                                                        "Gamma_male_young", "Gamma_female_young", 
                                                        "Mix_male_young",  "Mix_female_young", 
                                                        "Gamma_male_old", "Gamma_female_old", 
                                                        "Mix_male_old",  "Mix_female_old"))


## Low count filtering
left_ventrical_raw_counts <- left_ventrical_raw_counts[,lv_coldata$sample_id]
left_ventrical_raw_counts_filtered <- left_ventrical_raw_counts[which(rowSums(left_ventrical_raw_counts >= 10) >= 50),]



## LV female DE analysis
lv_rad_coldata_females <- lv_coldata[which(as.character(lv_coldata$gender) == "female"),]

lv_rad_coldata_females$group <- factor(lv_rad_coldata_females$group, levels = c("Sham_female_young", "Sham_female_old", "Mix_female_young", "Mix_female_old", "Gamma_female_young", "Gamma_female_old"
))

lv_rad_female_counts <- left_ventrical_raw_counts_filtered[,lv_rad_coldata_females$sample_id]


lv_dds_females <- DESeqDataSetFromMatrix(countData = lv_rad_female_counts, 
                                         colData = lv_rad_coldata_females, 
                                         design = ~Field+time_collect)
lv_dds_females <- DESeq(lv_dds_females)


## LV male DE analysis
lv_rad_coldata_males <- lv_coldata[which(as.character(lv_coldata$gender) == "male"),]
lv_rad_coldata_males$group <- factor(lv_rad_coldata_males$group, levels = c("Sham_male_young", "Sham_male_old", "Mix_male_young", "Mix_male_old", "Gamma_male_young", "Gamma_male_old"))

lv_rad_male_counts <- left_ventrical_raw_counts_filtered[,lv_rad_coldata_males$sample_id]


lv_dds_males <- DESeqDataSetFromMatrix(countData = lv_rad_male_counts, 
                                       colData = lv_rad_coldata_males, 
                                       design = ~Field+time_collect)
lv_dds_males <- DESeq(lv_dds_males)



lv_rad_deg_tables <- c(
  lapply(list(Female_Gamma = as.data.frame(results(lv_dds_females, contrast=c("Field","Gamma","Sham"))),
              Female_Mix = as.data.frame(results(lv_dds_females, contrast=c("Field","Mix","Sham")))), function(x) {
                x$gene_name <- gene_id_to_symbol[rownames(x)]
                x$gene_type <- gene_id_to_type[rownames(x)]
                x$gene_description <- gene_id_to_description[rownames(x)]
                x$chr <- gene_id_to_chromosome[rownames(x)]
                x[order(x$padj),]
              }),
  lapply(list(Male_Gamma = as.data.frame(results(lv_dds_males, contrast=c("Field","Gamma","Sham"))),
              Male_Mix = as.data.frame(results(lv_dds_males, contrast=c("Field","Mix","Sham")))), function(x) {
                x$gene_name <- gene_id_to_symbol[rownames(x)]
                x$gene_type <- gene_id_to_type[rownames(x)]
                x$gene_description <- gene_id_to_description[rownames(x)]
                x$chr <- gene_id_to_chromosome[rownames(x)]
                x[order(x$padj),]
              })
)




### Upset plots of left and right ventricle DEG overlaps
## Feamles
upset(
  fromList(setNames(nm = c("LV_Female_Gamma", "LV_Female_simGCRsim", "RV_Female_Gamma", "RV_Female_simGCRsim"), object = c(
    setNames(nm = paste0("LV_",names(lv_rad_deg_tables)), 
             object = lapply(lv_rad_deg_tables, function(x) {
               rownames(x)[which(x$padj < 0.05)]
             })
    )[1:2],
    setNames(nm = paste0("RV_",names(rv_rad_deg_tables)), 
             object = lapply(rv_rad_deg_tables, function(x) {
               rownames(x)[which(x$padj < 0.05)]
             })
    )[1:2]
  ))), 
  keep.order = T,
  # order.by = "freq", 
  matrix.color = "#228b22", set_size.show = TRUE,
  point.size = 5, line.size = 1.5, nintersects = NA, text.scale = c(2, 2, 2, 2, 2.5, 3), nsets = 4, set_size.scale_max = 45, mainbar.y.max = 400 
)


## Males
upset(
  fromList(setNames(nm = c("LV_Male_Gamma", "LV_Male_simGCRsim", "RV_Male_Gamma", "RV_Male_simGCRsim"), object = c(
    setNames(nm = paste0("LV_",names(lv_rad_deg_tables)), 
             object = lapply(lv_rad_deg_tables, function(x) {
               rownames(x)[which(x$padj < 0.05)]
             })
    )[3:4],
    setNames(nm = paste0("RV_",names(rv_rad_deg_tables)), 
             object = lapply(rv_rad_deg_tables, function(x) {
               rownames(x)[which(x$padj < 0.05)]
             })
    )[3:4]
  ))), 
  keep.order = T,
  # order.by = "freq", 
  matrix.color = "#228b22", set_size.show = TRUE,
  point.size = 5, line.size = 1.5, nintersects = NA, text.scale = c(2, 2, 2, 2, 2.5, 3), nsets = 4, set_size.scale_max = 450
)
