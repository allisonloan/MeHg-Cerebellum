#Data Import and Filtration---------

#Import Data
library(readxl)
library(dplyr)
Batch1_ori <- read_excel("Tuana_proteinGroups.xlsx")

#Remove any proteins with + in potential contaminants (only values in column are NA and +)

Batch1_ori$`Potential contaminant` <- as.character(Batch1_ori$`Potential contaminant`)

Batch1_ori <- Batch1_ori %>%
  filter(is.na(`Potential contaminant`))

unique(Batch1_ori$`Potential contaminant`)

#light filter 

Batch1_LFQ<-Batch1_ori[,c(1,68:75)]

newnames<- c("Protein IDs", "Control_Male_1", "Control_Male_2", 
             "Control_Female_1", "Control_Female_2", 
             "Mehg_Male_1", "Mehg_Male_2",
             "Mehg_Female_1", "Mehg_Female_2")

colnames(Batch1_LFQ) <- newnames
rm(newnames,Batch1_ori)


Batch1_LFQ2<- Batch1_LFQ
Batch1_LFQ2[Batch1_LFQ2 == 0] <- NA

Na_count <- rowSums(is.na(Batch1_LFQ2[, -1]))  # Count NAs across all samples (excluding protein ID)

# Filter out proteins with NAs in at least 6/8 samples
Mouse3_LFQ<- Batch1_LFQ2[Na_count < 6, ]

rm(Na_count,Batch1_LFQ)

library(stringr)
Mouse3_LFQ <- Mouse3_LFQ %>%
  mutate(`Protein IDs` = sapply(str_split(`Protein IDs`, ";"), `[`, 1))


#Venn Diagram-------

library(ggvenn)

Control <- c("Control_Male_1", "Control_Male_2", "Control_Female_1", "Control_Female_2")  
Mehg <- c("Mehg_Male_1", "Mehg_Male_2", "Mehg_Female_1", "Mehg_Female_2") 

all_Control <- Mouse3_LFQ %>%
  filter(
    rowSums(!is.na(select(., all_of(Control)))) >= 2,  # Present in at least 2/4 Control
  )

all_Mehg <- Mouse3_LFQ %>%
  filter(
    rowSums(!is.na(select(., all_of(Mehg)))) >= 2,  # Present in at least 2/4 Mehg
  )


venn_data<- list(
  Control = all_Control$`Protein IDs`,
  Mehg = all_Mehg$`Protein IDs`)

ggvenn(venn_data, c("Control", "Mehg"), fill_color=c('#FFFFC5', "#89CFF0"))


# Find unique proteins for each group
unique_Control <- setdiff(all_Control$`Protein IDs`, all_Mehg$`Protein IDs`)


unique_Mehg <- setdiff(all_Mehg$`Protein IDs`, all_Control$`Protein IDs`)

# Find common proteins
common_proteins <- intersect(all_Control$`Protein IDs`, all_Mehg$`Protein IDs`)

names(common) <- gsub("sp\\|([A-Z0-9]+)\\|.*", "\\1", names(common))

Unique_Control<- Mouse3_LFQ[unique_Control,]
Unique_Mehg<- Mouse3_LFQ[unique_Mehg,]

Unique_Mehg <- Mouse_clean %>%
  filter(`Protein IDs` %in% unique_Mehg)

Unique_Mehg <- data.frame(Unique_Mehg, row.names = 1)

rownames(Unique_Mehg) <- gsub("sp\\|([A-Z0-9]+)\\|.*", "\\1", rownames(Unique_Mehg))

selected_ids <- c(unique_Control, unique_Mehg, common_proteins)

#Mouse clean contains all unique and common proteins after filtering
Mouse_clean <- Mouse3_LFQ %>%
  filter(`Protein IDs` %in% selected_ids)

Common_venn <- Mouse_clean %>%
  filter(`Protein IDs` %in% common_proteins)

write.xlsx(Unique_Mehg, "Unique_Mehg.xlsx")

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "UniqueControl")
writeData(wb, sheet = "UniqueControl", Unique_Control,rowNames = TRUE)
addWorksheet(wb, "UniqueMehg")
writeData(wb, sheet = "UniqueMehg", Unique_Mehg, rowNames = TRUE)
addWorksheet(wb, "Common")
writeData(wb, sheet="Common", Common_venn, rowNames = TRUE)
saveWorkbook(wb, file = "Venn Output.xlsx", overwrite = TRUE)

rm(wb, unique_Control,unique_Mehg,common_proteins,venn_data,all_Control, selected_ids,
   Control,Mehg,Unique_Mehg,Unique_Control,all_Mehg,Batch1_LFQ2)


#Data Visualization----
# Plot all protein intensities of one sample to see distribution

Mouse_clean <- data.frame(Mouse_clean, row.names = 1)

b2_exc <- b2_exc %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()

b2_exc<-b2_exc[-1,]

b2_exc[] <- lapply(b2_exc, as.numeric)

# Plot all protein intensities of one sample
plot(x = Mouse_clean$`Control_Male_2`)

# Log10 of variable intensity
plot(x = log2(Mouse_clean$`Control_Male_2`))

# Plot all protein intensities of one sample as histogram
hist(x = Mouse_clean$`Control_Male_1`)

# Log10 of variable intensity
hist(x = log2(Mouse_clean$`Mehg_Female_1`), 20)


# More Data Visualization ---- 

#visualizing the data 

# Compare Biological Replicates
plot(x = Batch2_filtered$`LFQ intensity Mehg_Male_2A`, 
     y = Batch2_filtered$`LFQ intensity Mehg_Male_3A`)

# Log10
plot(x = log2(Mouse_clean$`Mehg_Male_1`), 
     y = log2(Mouse_clean$`Mehg_Male_3`))


#PCA--------

Mouse_imp <- data.frame(Mouse_clean, row.names = 1)

#proteins as columns, samples as rows
Mouse_imp <- Mouse_imp %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()


Mouse_imp.log2 <- log2(Mouse_imp)

Mouse_imp[is.na(Mouse_imp)] <- 0

shift <- 1.8
width <- 0.3
# for loop iterates through all columns that contain missing values

set.seed(100)

for (i in which(apply(X = Mouse_imp, MARGIN = 2, FUN = function(x) any(x == 0)))) {
  
  # Impute from normal distribution rnorm (based on log2 values)
  Mouse_imp[Mouse_imp[, i] == 0, i] <- round(2 ^ rnorm(n = sum(Mouse_imp[, i] == 0),
                                                       mean = mean(Mouse_imp.log2[, i], na.rm = TRUE) - 
                                                         shift * sd(Mouse_imp.log2[, i], na.rm = TRUE),
                                                       sd = width * sd(Mouse_imp.log2[, i], na.rm = TRUE)),
                                             digits = -2)
  
}

rm(Mouse_imp.log2, i,shift,width)
log2.Mouse_imp <- log2(Mouse_imp)

# Center and Scale data
pca_scaled <- log2.Mouse_imp %>%
  dplyr::mutate(dplyr::across(.cols = where(is.numeric),
                              .fns = function(x) c(scale(x,scale=TRUE))))

pca_scaled <- log2.Mouse_imp %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ scale(.)))

# Compute principal components with prcomp() function
data_prcomp <- pca_scaled  %>%
  dplyr::select(where(is.numeric)) %>%
  prcomp()

gender <- factor(x = c("male", "male","female", "female",
                       "male","male","female","female"), 
                 levels = c("male", "female"))

groups <- factor(x = c("control", "control", "control","control",
                       "mehg","mehg","mehg", "mehg"), 
                 levels = c("control", "mehg"))

# Add observations names to results data frame
data_pca <- pca_scaled %>%
  dplyr::select(-where(is.numeric)) %>% 
  dplyr::mutate(Groups = groups) %>%
  cbind(data_prcomp[["x"]]) %>% 
  as_tibble()


# Plot PCA

x <- "PC1"
y <- "PC2"


p <- ggplot(data_pca, aes(x = .data[[x]],
                          y = .data[[y]],
                          color = groups,
                          fill = groups)) +
  geom_point() +
  theme_minimal()+
  coord_cartesian(ylim=c(-50,50), xlim=c(-70,70))
#geom_hline(yintercept=c(0),col="black", linetype="dashed")+
#geom_vline(xintercept=c(0),col="black", linetype="dashed")

# Add fraction of variance
p <- p +
  xlab(paste0(
    x,
    " (",
    round(
      100 *
        data_prcomp[["sdev"]][as.numeric(substring(x, 3))]^2 /
        sum(data_prcomp[["sdev"]]^2),
      digits = 1), "%)")) +
  ylab(paste0(
    y,
    " (",
    round(
      100 *
        data_prcomp[["sdev"]][as.numeric(substring(y, 3))]^2 /
        sum(data_prcomp[["sdev"]]^2),
      digits = 1), "%)"))

p



#environment cleanup
rm(pca_scaled,p,data_prcomp,data_pca, i, x,y)


#2Way ANOVA----

log2.common <- log2.common %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()

Anova_data<-log2.common


Anova_data2<-setDT(Anova_data, keep.rownames = "Observations")[]

Anova_data2 <- Anova_data2 %>%
  mutate(Sex = gender) %>%
  select(Observations,Sex, everything())

Anova_data2 <- Anova_data2 %>%
  mutate(Groups = groups) %>%
  select(Observations,Groups,Sex, everything())

# List to save models in
list_aov <- list()

# Define formula
formula <- "x ~ Groups*Sex"

# Formulate expression for model
expr.left <- substr(formula, 1, regexpr("x", formula) - 1)
expr.right <- substring(formula, regexpr("x", formula) + 1)

names(Anova_data2) <- gsub("sp\\|([A-Z0-9]+)\\|.*", "\\1", names(Anova_data2))


# Apply aov model to variables
for (i in colnames(Anova_data2)[-c(1, 2,3)]) {
  # Combine protein name to formula
  expr <- paste0(expr.left, i, expr.right)
  # Do ANOVA 
  list_aov[[i]] <- aov(formula = rlang::eval_tidy(rlang::parse_expr(expr)), 
                       data = Anova_data2)
  
}

# Extract p-values
list_aov_p <- lapply(list_aov, function(x) {
  summary_x <- summary(x)
  
  # Extract all p-values
  p_values <- summary_x[[1]][, "Pr(>F)"]
  
  # Assign names based on row names
  names(p_values) <- rownames(summary_x[[1]])
  
  # Exclude residuals
  p_values <- p_values[!names(p_values) %in% "Residuals"]
  
  return(p_values)
})

p_values_df <- do.call(rbind, list_aov_p)
p_values_df<-p_values_df[,-4]

# Post-hoc test
list_posthoc <- lapply(list_aov, FUN = TukeyHSD)

# Extract TukeyHSD p-values
list_posthoc_p <- lapply(list_posthoc, function(x) {
  lapply(x, function(comp) {
    p_adj_values <- comp[, "p adj"]
    data.frame(t(p_adj_values))
  })
})

extract_interaction_p_values <- function(p_values_list) {
  # Assuming the first two elements are Factor 1 and Factor 2 p-values,
  # and the remaining 6 elements are the interaction p-values.
  interaction_p_values <- p_values_list[3:8]
  
  # Convert to a named vector
  return(unlist(interaction_p_values))
}

# Apply this function to each element (variable) in the list
interaction_p_values_list <- lapply(list_posthoc_p, extract_interaction_p_values)

# Combine the results into a data frame
interaction_p_values_df <- do.call(rbind, interaction_p_values_list)

total_anova_turkey_pvalues<-cbind(p_values_df,interaction_p_values_df)

#extract the diff values from posthoc
list_posthoc_diff <- lapply(list_posthoc, function(x) {
  lapply(x, function(comp) {
    diff_values <- comp[, "diff"]
    data.frame(t(diff_values))
  })
})

extract_interaction_diff_values <- function(diff_values_list) {
  # Assuming the first two elements are Factor 1 and Factor 2 p-values,
  # and the remaining 6 elements are the interaction p-values.
  interaction_diff_values <- diff_values_list[3:8]
  
  # Convert to a named vector
  return(unlist(interaction_diff_values))
}

# Apply this function to each element (variable) in the list
interaction_diff_values_list <- lapply(list_posthoc_diff, extract_interaction_diff_values)


# Combine the results into a data frame
interaction_diff_values_df <- do.call(rbind, interaction_diff_values_list)

colnames(interaction_diff_values_df) <- gsub(pattern = "Groups:Sex", 
                                             replacement = "log2fc_", 
                                             x = colnames(interaction_diff_values_df))

# Combine p-values with fold-changes
data_anova_results <- cbind(total_anova_turkey_pvalues,interaction_diff_values_df)


# Output data
data_anova_results <- data_anova_results %>% 
  dplyr::mutate(`Protein names` = data_prot$`Protein names`, 
                `Gene names` = data_prot$`Gene names`, 
                .after = variables)

# Export table
openxlsx::write.xlsx(data_anova_results, "ANOVA_results.xlsx",rowNames=TRUE)

#Clean Environment
rm(total_anova_turkey_pvalues,p_values_df,list_posthoc_p,list_posthoc_diff,list_posthoc,list_aov_p,list_aov,
   interaction_p_values_list, interaction_p_values_df,interaction_diff_values_list,interaction_diff_values_df)
rm(x,y,pAdjustMethod,p.value.cutoff,i,formula,expr.right, expr.left, expr)

#Volcano Plot After ANOVA------

# Volcano plot Control v Mehg Male
controlvmehg_male<-data_anova_results[,c(4,10)]
colnames(controlvmehg_male) <- c("p.value", "log2fc")
#just pvalue threshold
controlvmehg_male <- as.data.frame(controlvmehg_male)
controlvmehg_male$diffexp<-'NO'
controlvmehg_male$diffexp[controlvmehg_male$log2fc> log2(1.5) & controlvmehg_male$p.value<0.05]<-'UP'
controlvmehg_male$diffexp[controlvmehg_male$log2fc<  log2(1/1.5) & controlvmehg_male$p.value<0.05]<-'DOWN'


ggplot(data = controlvmehg_male,
       aes(x = log2fc,
           y = -log10(p.value),
           col=diffexp)) +
  #geom_hline(yintercept=c(-log10(0.05)),col="black", linetype="dashed")+
  #geom_vline(xintercept=c(-0.6,0.6),col="black", linetype="dashed")+
  coord_cartesian(ylim=c(-0,4.1), xlim=c(-4,4))+
  geom_point(shape = 16, stroke = 0, size=2) +
  theme_classic() +
  scale_color_manual(values = c("UP" = "red",
                                "DOWN" = "blue",
                                "NO" = "grey"),
                     labels=c("Downregulated","Not Significant", "Up Regulated")) +
  labs(
    title = "Volcano Plot of Control Male v Mehg Male",
    x = "Log2 Fold Change",
    y = "-Log10(p-value)",
    color = "Significance"
  ) 



# Volcano plot Control Female v Mehg Female
controlvmehg_female<-data_anova_results[,c(9,15)]
colnames(controlvmehg_female) <- c("p.value", "log2fc")
controlvmehg_female <- as.data.frame(controlvmehg_female)
controlvmehg_female$diffexp<-'NO'
controlvmehg_female$diffexp[controlvmehg_female$log2fc> log2(1.5) & controlvmehg_female$p.value<0.05]<-'UP'
controlvmehg_female$diffexp[controlvmehg_female$log2fc< log2(1/1.5) &controlvmehg_female$p.value<0.05]<-'DOWN'


ggplot(data = controlvmehg_female,
       aes(x = log2fc,
           y = -log10(p.value),
           col=diffexp)) +
  #geom_hline(yintercept=c(-log10(0.05)),col="black", linetype="dashed")+
  #geom_vline(xintercept=c(-0.6,0.6),col="black", linetype="dashed")+
  coord_cartesian(ylim=c(0,4), xlim=c(-6,6))+
  geom_point(shape = 16, stroke = 0, size=2) +
  theme_classic() +
  scale_color_manual(values = c("UP" = "red",
                                "DOWN" = "blue",
                                "NO" = "grey"),
                     labels=c("Downregulated","Not Significant", "Up Regulated")) +
  labs(
    title = "Volcano Plot of Control Female v Mehg Female",
    x = "Log2 Fold Change",
    y = "-Log10(p-value)",
    color = "Significance"
  ) 

anova_results_male <- data.frame(ProteinID = rownames(controlvmehg_male), 
                                 controlvmehg_male, row.names = NULL)


anova_results_female <- data.frame(ProteinID = rownames(controlvmehg_female), 
                                   controlvmehg_female, row.names = NULL)


combined_volcano <- rbind(
  data.frame(ProteinID = anova_results_male$ProteinID,
             log2FC = anova_results_male$log2fc,
             p.value = anova_results_male$p.value,
             group = "Male",
             diffexp = anova_results_male$diffexp),
  
  data.frame(ProteinID = anova_results_female$ProteinID,
             log2FC = anova_results_female$log2fc,
             p.value = anova_results_female$p.value,
             group = "Female",
             diffexp = anova_results_female$diffexp)
)

ggplot(data = combined_volcano,
       aes(x = log2FC,
           y = -log10(p.value),
           col=diffexp,
           shape = group)) +
  #geom_hline(yintercept=c(-log10(0.05)),col="black", linetype="dashed")+
  #geom_vline(xintercept=c(-0.6,0.6),col="black", linetype="dashed")+
  coord_cartesian(ylim=c(-0,4), xlim=c(-4,4))+
  geom_point(stroke = 0, size=2) +
  theme_classic() +
  scale_color_manual(values = c("UP" = "red",
                                "DOWN" = "blue",
                                "NO" = "grey"),
                     labels=c("Downregulated","Not Significant", "Up Regulated")) +
  labs(
    title = "Volcano Plot of Control vs. Mehg ",
    x = "Log2 Fold Change",
    y = "-Log10(p adj-value)",
    color = "Significance"
  ) 

#Export Volcano Plot Data
wb <- createWorkbook()
addWorksheet(wb, "Control Female-Mehg Female")
writeData(wb, sheet = "Control Female-Mehg Female", controlvmehg_female, rowNames = TRUE)
addWorksheet(wb, "Control Male-Mehg Male")
writeData(wb, sheet = "Control Male-Mehg Male", controlvmehg_male, rowNames = TRUE)
saveWorkbook(wb, file = "VolcanoPlots_data.xlsx", overwrite = TRUE)

#Diff Exp Heatmap Final(After ANOVA)----

#results_df2 <- data.frame(results_df, row.names = 1)

heatmap.c <- list(controlvmehg_female,controlvmehg_male)

diffexp_ids <- unique(unlist(lapply(heatmap.c, function(df) {
  rownames(df)[df$diffexp%in% c("UP","DOWN")]
})))

names(log2.Mouse_imp) <- gsub("sp\\|([A-Z0-9]+)\\|.*", "\\1", names(log2.Mouse_imp))

heatmap_df<-log2.Mouse_imp[,diffexp_ids]

hmap_data_scaled <- heatmap_df %>%
  dplyr::mutate(dplyr::across(.cols = where(is.numeric),
                              .fns = function(x) c(scale(x))))

# Transform data frame into matrix format
hmap_data_scaled <- t(hmap_data_scaled)

library(pheatmap)
pheatmap::pheatmap(mat = hmap_data_scaled,show_rownames=TRUE,
                   clustering_distance_rows="euclidean",cluster_cols = TRUE,
                   fontsize_row = 10,
                   main="Diff. Expression; Up and Down Regulated Proteins in Mehg Exposed Mice")

hmap_all<- log2.Mouse_imp[,diffexp_ids] #diffexp_ids is a vector containing all diff exp ids of both m and f after anova
hmap_male<- log2.Mouse_imp[,mdiffexp_ids]# mdiffexp_ids only male diffexp ids after anova
hmap_female<- log2.Mouse_imp[,fdiffexp_ids] #fdiffexp_ids only female diffexp ids after anova

wb <- createWorkbook()
addWorksheet(wb, "HeatMap_Male")
writeData(wb, sheet = "HeatMap_Male", hmap_male, rowNames = TRUE)
addWorksheet(wb, "HeatMap_Female")
writeData(wb, sheet = "HeatMap_Female", hmap_female, rowNames = TRUE)
addWorksheet(wb, "HeatMap_All")
writeData(wb, sheet = "HeatMap_All", hmap_all, rowNames = TRUE)
saveWorkbook(wb, file = "HmapData.xlsx", overwrite = TRUE)

rm(wb,i,contrast_matrix,design,data_t.test,fit,fit2,list_t.test,metadata,results)

#Kegg(After ANOVA)----

kegg.c <- list(controlvmehg_female)

kegg.c <- list(controlvmehg_male,controlvmehg_female)

diffexp_ids <- unique(unlist(lapply(kegg.c, function(df) {
  rownames(df)[df$diffexp %in% c("UP","DOWN")]
})))

unique_Mehg<- rownames(Unique_Mehg)



unique_m_mehg<-c("A2A432","O35526","P47758","Q61282", "Q68ED7","Q80X71","Q8BFR4")
selected_ids <- c(diffexp_ids,unique_Mehg)

library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
id_conversion <- getBM(
  attributes = c("uniprot_gn_id", "entrezgene_id"),
  filters = "uniprot_gn_id",
  values = selected_ids,
  mart = mart
)

entrez_ids<-id_conversion$entrezgene_id

#perform the enrichment  

kegg_enrich <- enrichKEGG(gene = entrez_ids, 
                          organism = "mmu",     
                          pAdjustMethod = "none",
                          qvalueCutoff = 0.05)

kegg_tog<- as.data.frame(kegg_enrich@result)

kegg_tog <- kegg_tog %>% mutate(
  geneRatioDecimal = sapply(GeneRatio, function(x) {
    nums <- as.numeric(unlist(strsplit(x, "/")))
    return(nums[1] / nums[2])
  })
)

# Print the result to check
print(kegg_all)

ggplot(kegg_combined, aes(x = geneRatioDecimal, y = reorder(Description, geneRatioDecimal), color = qvalue, size = Count)) +
  geom_point() +
  labs(title = "KEGG Enrichment Analysis", x = "Gene Ratio", y = "Term", color = "qvalue", size = "Gene Count") +
  theme_minimal()

openxlsx::write.xlsx(sig_pathways, "KEGG.xlsx",rowNames=TRUE)

wb <- createWorkbook()
addWorksheet(wb, "KEGG_diffexp BH")
writeData(wb, sheet = "KEGG_diffexp BH", kegg_diffexp,rowNames = TRUE)
addWorksheet(wb, "KEGG_diffexp None")
writeData(wb, sheet = "KEGG_diffexp None", kegg_diffexp2,rowNames = TRUE)
addWorksheet(wb, "KEGG Mehg Unique None")
writeData(wb, sheet = "KEGG Mehg Unique None", kegg_unique, rowNames = TRUE)
addWorksheet(wb, "Diff+Unique_RunTog")
writeData(wb, sheet = "Diff+Unique_RunTog", kegg_tog, rowNames = TRUE)
addWorksheet(wb, "All Sig Pathways")
writeData(wb, sheet = "All Sig Pathways", kegg_combined, rowNames = TRUE)
saveWorkbook(wb, file = "AllKEGGFinal.xlsx", overwrite = TRUE)

#GO (After ANOVA)----
library(GOstats)

selected_ids <- c(unique_Mehg,diffexp_ids)

library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
id_conversion <- getBM(
  attributes = c("uniprot_gn_id", "entrezgene_id"),
  filters = "uniprot_gn_id",
  values = selected_ids,
  mart = mart
)

#make vector with new entrez ids
entrez_ids<-id_conversion$entrezgene_id

go_BP2 <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)

bp_df <- as.data.frame(go_BP)
bp_df2 <- as.data.frame(go_BP2)

BP_upreg <- enrichGO(
  gene = upreg_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
bp_upreg <- as.data.frame(BP_upreg)

BP_downreg <- enrichGO(
  gene = downreg_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
bp_downreg <- as.data.frame(BP_downreg)

BP_unique <- enrichGO(
  gene = unique_Mehg,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)

bp_unique <- as.data.frame(BP_unique)


go_CC2 <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)

cc_df <- as.data.frame(go_CC)
cc_df2 <- as.data.frame(go_CC2)

CC_upreg <- enrichGO(
  gene = upreg_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
cc_upreg <- as.data.frame(CC_upreg)

CC_downreg <- enrichGO(
  gene = downreg_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
cc_downreg <- as.data.frame(CC_downreg)

CC_unique <- enrichGO(
  gene = unique_Mehg,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)

cc_unique <- as.data.frame(CC_unique)

go_MF2 <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
mg_df2 <- as.data.frame(go_MF2)

MF_upreg <- enrichGO(
  gene = upreg_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
mf_upreg <- as.data.frame(MF_upreg)

MF_downreg <- enrichGO(
  gene = downreg_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
mf_downreg <- as.data.frame(MF_downreg)

MF_unique <- enrichGO(
  gene = unique_Mehg,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Use "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
  pAdjustMethod = "none",  # p-value adjustment method
  qvalueCutoff = 0.05  # Adjust threshold for significance
)
MF_unique <- as.data.frame(MF_unique)


ggplot(bp_df2[1:20,], aes(x = reorder(Description, Count), y = Count, fill=qvalue)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "20 BP GO Enrichment", x = "Term", y = "Gene Count") +
  theme_minimal() +
  coord_flip() 

ggplot(cc_df2[1:20,], aes(x = reorder(Description, Count), y = Count, fill=qvalue)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "20 CC GO Enrichment", x = "Term", y = "Gene Count") +
  theme_minimal() +
  coord_flip() 

ggplot(mg_df2[1:20,], aes(x = reorder(Description, Count), y = Count, fill=qvalue)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "20 MF GO Enrichment", x = "Term", y = "Gene Count") +
  theme_minimal() +
  coord_flip() 

# Get the top 5 terms for each category
top_BP <- head(BP_upreg, 5)
top_BP$ONTOLOGY<- "BP"
top_CC <- head(cc_upreg, 5)
top_CC$ONTOLOGY<- "CC"
top_MF <- head(MF_upreg, 5)
top_MF$ONTOLOGY<- "MF"

# Combine the top 5 terms from all categories
top_GO_upreg <- rbind(top_BP, top_CC, top_MF)
top_GO_upreg <- top_GO_upreg %>%
  mutate(Description = factor(Description, levels = unique(Description[order(ONTOLOGY, qvalue)])))

ggplot(top_GO_upreg, aes(x = Description, y = -log10(qvalue), fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(aes(label = Count), hjust = 1.5, size = 14) +  # Display count above bars
  labs(title = "Top 5 GO Terms Upregulated by Category (BP, CC, MF)",
       x = "GO Term", y = "-log10(p.adjust)", fill = "Ontology") +
  theme_minimal(base_family = "Arial", base_size = 30) +
  coord_flip() +  # Flip axes for better readability
  scale_fill_manual(values = c("BP" = "#f0f921", "MF" = "#f89540", "CC" = "#cc4778")) +  # Custom colors for ONTOLOGY
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.7),
    legend.position = "right",
    axis.text.y = element_text(size=25),
    axis.text.x = element_text(size=25)
  )

#exported bar plot as width 900, height 600

# Get the top 5 terms for each category
top_BP <- head(BP_downreg, 5)
top_BP$ONTOLOGY<- "BP"
top_CC <- head(cc_downreg, 5)
top_CC$ONTOLOGY<- "CC"
top_MF <- head(MF_downreg, 5)
top_MF$ONTOLOGY<- "MF"

# Combine the top 5 terms from all categories
top_GO_downreg <- rbind(top_BP, top_CC, top_MF)
top_GO_downreg <- top_GO_downreg %>%
  mutate(Description = factor(Description, levels = unique(Description[order(ONTOLOGY, qvalue)])))

ggplot(top_GO_downreg, aes(x = Description, y = -log10(qvalue), fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(aes(label = Count), hjust = 1.5, size = 14) +  # Display count above bars
  labs(title = "Top 5 GO Terms Downregulated by Category (BP, CC, MF)",
       x = "GO Term", y = "-log10(p.adjust)", fill = "Ontology") +
  theme_minimal(base_family = "Arial", base_size = 30) +
  coord_flip() +  # Flip axes for better readability
  scale_fill_manual(values = c("BP" = "#f0f921", "MF" = "#f89540", "CC" = "#cc4778")) +  # Custom colors for ONTOLOGY
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.7),
    legend.position = "right",
    axis.text.y = element_text(size=25),
    axis.text.x = element_text(size=25)
  )
