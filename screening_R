#############################################################
#### Script for the whole analysis of screening results #####
#############################################################

## upload libraries
library(tidyverse)
library(dplyr)
library(readxl)
library(writexl)
library(reshape2)
library(MESS)
library(ggplot2)
library(scales)
library(svglite)
library(GO.db)

########################################################################

#read the most recent excel
#c_temp <- read_excel("./screening_data_with_genes.xlsx")

#add columns with measurement value minus T0 measurement value
c_temp$gene_CTR_T0_minusT0 <- c_temp$gene_CTR_T0 - c_temp$gene_CTR_T0
c_temp$gene_CTR_T1_minusT0 <- c_temp$gene_CTR_T1 - c_temp$gene_CTR_T0
c_temp$gene_CTR_T2_minusT0 <- c_temp$gene_CTR_T2 - c_temp$gene_CTR_T0
c_temp$gene_CTR_T3_minusT0 <- c_temp$gene_CTR_T3 - c_temp$gene_CTR_T0
c_temp$gene_CTR_T4_minusT0 <- c_temp$gene_CTR_T4 - c_temp$gene_CTR_T0
c_temp$gene_TAT_T0_minusT0 <- c_temp$gene_TAT_T0 - c_temp$gene_TAT_T0
c_temp$gene_TAT_T1_minusT0 <- c_temp$gene_TAT_T1 - c_temp$gene_TAT_T0
c_temp$gene_TAT_T2_minusT0 <- c_temp$gene_TAT_T2 - c_temp$gene_TAT_T0
c_temp$gene_TAT_T3_minusT0 <- c_temp$gene_TAT_T3 - c_temp$gene_TAT_T0
c_temp$gene_TAT_T4_minusT0 <- c_temp$gene_TAT_T4 - c_temp$gene_TAT_T0

c_temp$WT_CTR_T0_minusT0 <- c_temp$WT_CTR_T0 - c_temp$WT_CTR_T0
c_temp$WT_CTR_T1_minusT0 <- c_temp$WT_CTR_T1 - c_temp$WT_CTR_T0
c_temp$WT_CTR_T2_minusT0 <- c_temp$WT_CTR_T2 - c_temp$WT_CTR_T0
c_temp$WT_CTR_T3_minusT0 <- c_temp$WT_CTR_T3 - c_temp$WT_CTR_T0
c_temp$WT_CTR_T4_minusT0 <- c_temp$WT_CTR_T4 - c_temp$WT_CTR_T0
c_temp$WT_TAT_T0_minusT0 <- c_temp$WT_TAT_T0 - c_temp$WT_TAT_T0
c_temp$WT_TAT_T1_minusT0 <- c_temp$WT_TAT_T1 - c_temp$WT_TAT_T0
c_temp$WT_TAT_T2_minusT0 <- c_temp$WT_TAT_T2 - c_temp$WT_TAT_T0
c_temp$WT_TAT_T3_minusT0 <- c_temp$WT_TAT_T3 - c_temp$WT_TAT_T0
c_temp$WT_TAT_T4_minusT0 <- c_temp$WT_TAT_T4 - c_temp$WT_TAT_T0

########## Fold change calculation ##########################
## Fold change gene growth without peptide compared to WT without peptide
c_temp$FoldChange_growth <- with(c_temp, gene_CTR_T4_minusT0 / mean(WT_CTR_T4_minusT0))

## Fold change TAT vs CTR at T3 and T4
c_temp$FoldChange_gene_T3 <- c_temp$gene_TAT_T3_minusT0 / c_temp$gene_CTR_T3_minusT0
c_temp$FoldChange_gene_T4 <- c_temp$gene_TAT_T4_minusT0 / c_temp$gene_CTR_T4_minusT0
c_temp$FoldChange_WT_T3 <- c_temp$WT_TAT_T3_minusT0 / c_temp$WT_CTR_T3_minusT0
c_temp$FoldChange_WT_T4 <- c_temp$WT_TAT_T4_minusT0 / c_temp$WT_CTR_T4_minusT0


#### Calculate foldchange gene vs WT
c_temp$FoldChange_genevsWT_T3 <- c_temp$FoldChange_gene_T3 / c_temp$FoldChange_WT_T3
c_temp$FoldChange_genevsWT_T4 <- c_temp$FoldChange_gene_T4 / c_temp$FoldChange_WT_T4

############################################################

############################################################
########### Setting thresholds for FoldChange and AUC ######

## calculate mean and sd

## T3
mean_FoldChange_T3_WT <- mean(c_temp$FoldChange_WT_T3)
sd_FoldChange_T3_WT <- sd(c_temp$FoldChange_WT_T3)
maxValue_T3 <- max(c_temp$FoldChange_gene_T3)

## T4
mean_FoldChange_T4_WT <- mean(c_temp$FoldChange_WT_T4)
sd_FoldChange_T4_WT <- sd(c_temp$FoldChange_WT_T4)

### plot the distribution of the points to select the best thresholds
c_temp_growth <- filter(c_temp, gene_CTR_T4_minusT0 > 0.4)

## T3 histogram
p <- ggplot(c_temp_growth, aes(FoldChange_gene_T3))
p + geom_histogram(fill = "lightblue", alpha = 0.5, size = 1, binwidth = 0.01) + 
  theme_minimal() + 
  geom_vline(xintercept = mean_FoldChange_T3_WT, color = "black", size = 1) +
  geom_text(aes(x = mean_FoldChange_T3_WT - 0.13, y = 4.8, label = "Mean WT"), color = "black", size = 5) +
  geom_text(aes(x = mean_FoldChange_T3_WT + 0.28, y = 4.8, label = "+ 2 * SD"), color = "black", size = 5) +
  geom_text(aes(x = 0.85, y = 4, label = "resistant"), color = "black", size = 5)+
  geom_vline(xintercept = mean_FoldChange_T3_WT + 2*sd_FoldChange_T3_WT, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Fold change (NG6h TAT-RasGAP / NG6h control)") +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
ggsave(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "/plots/T3_max_hist.png"), width = 8, height = 6)

## T3 density
p + annotate("rect", xmin = mean_FoldChange_T3_WT + 2*sd_FoldChange_T3_WT, xmax = maxValue_T3, ymin = 0, ymax = 5, alpha = 0.2) + 
  geom_density(fill = "lightblue", alpha = 0.5, size = 1) + 
  theme_minimal() + 
  geom_vline(xintercept = mean_FoldChange_T3_WT, color = "black", size = 1) +
  geom_text(aes(x = mean_FoldChange_T3_WT - 0.13, y = 4.8, label = "Mean WT"), color = "black", size = 5) +
  geom_text(aes(x = mean_FoldChange_T3_WT + 0.28, y = 4.8, label = "+ 2 * SD"), color = "black", size = 5) +
  geom_text(aes(x = 0.85, y = 4, label = "resistant"), color = "black", size = 5)+
  geom_vline(xintercept = mean_FoldChange_T3_WT + 2*sd_FoldChange_T3_WT, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Fold change (NG6h TAT-RasGAP / NG6h control)") +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
ggsave(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "/plots/T3_max_dens.png"), width = 8, height = 6)

## T3 density
p + geom_line(stat = "density", size = 1) + 
  theme_minimal() + 
  geom_vline(xintercept = mean_FoldChange_T3_WT, color = "black", size = 1) +
  geom_text(aes(x = mean_FoldChange_T3_WT - 0.13, y = 4.8, label = "Mean WT"), color = "black", size = 5) +
  geom_text(aes(x = mean_FoldChange_T3_WT + 0.28, y = 4.8, label = "+ 2 * SD"), color = "black", size = 5) +
  geom_text(aes(x = 0.85, y = 4, label = "resistant"), color = "black", size = 5)+
  geom_vline(xintercept = mean_FoldChange_T3_WT + 2*sd_FoldChange_T3_WT, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Fold change (NG6h TAT-RasGAP / NG6h control)") +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  scale_y_continuous(labels = percent, name = "percent")
ggsave(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "/plots/T3_max.png"), width = 8, height = 6)

## T3 density with arrow
p + geom_density(fill = "lightblue", alpha = 0.5, size = 1) + 
  theme_bw() + 
  geom_vline(xintercept = mean_FoldChange_T3_WT, color = "black", size = 1) +
  geom_text(aes(x = mean_FoldChange_T3_WT - 0.13, y = 4.8, label = "Mean WT"), color = "black", size = 5) +
  geom_text(aes(x = mean_FoldChange_T3_WT + 0.28, y = 4.8, label = "+ 2 * SD"), color = "black", size = 5) +
  geom_text(aes(x = 0.85, y = 1.2, label = "resistant"), color = "black", size = 5)+
  geom_vline(xintercept = mean_FoldChange_T3_WT + 2*sd_FoldChange_T3_WT, linetype = "dashed", color = "black", size = 1) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13)) +
  geom_segment(aes(x = mean_FoldChange_T3_WT + 2*sd_FoldChange_T3_WT, y = 1, xend = 1.0, yend = 1), arrow = arrow(length = unit(0.5, "cm")), size = 1) +
  xlab("Fold change (NG6h TAT-RasGAP317-326 / NG6h control)")
ggsave(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "/plots/T3_arrow.svg"), width = 6, height = 4)

## T4
q <- ggplot(c_temp_growth, aes(FoldChange_gene_T4))
q + geom_density(fill = "red", alpha = 0.5, size = 1) + 
  theme_bw() + 
  geom_vline(xintercept = mean_FoldChange_T4_WT, color = "black", size = 1) +
  geom_text(aes(x = mean_FoldChange_T4_WT + 0.16, y = 3.1, label = "Mean WT"), color = "black", size = 5) +
  geom_vline(xintercept = mean_FoldChange_T4_WT - 3*sd_FoldChange_T4_WT, linetype = "dashed", color = "black", size = 1) +
  geom_text(aes(x = mean_FoldChange_T4_WT - 0.32, y = 3.1, label = "- 1*SD"), color = "black", size = 5) +
  geom_vline(xintercept = mean_FoldChange_T4_WT - sd_FoldChange_T4_WT, linetype = "dashed", color = "black", size = 1) +
  geom_text(aes(x = mean_FoldChange_T4_WT - 0.70, y = 3.1, label = "- 3*SD"), color = "black", size = 5)+
  labs(x = "Fold Change (NG24h TAT-RasGAP317-326 / NG24h control)") +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13)) +
  geom_segment(aes(x = mean_FoldChange_T4_WT - 3*sd_FoldChange_T3_WT, y = 1, xend = -0.15, yend = 1), arrow = arrow(length = unit(0.4, "cm")), size = 1) +
  geom_text(aes(x = -0.01, y = 1.2, label = "hypersensitive"), color = "black", size = 5) +
  scale_x_continuous(limits = c(-0.15, 1.5)) +
  geom_text(aes(x = mean_FoldChange_T4_WT - 0.39, y = 2.4, label = "more"), color = "black", size = 5) +
  geom_text(aes(x = mean_FoldChange_T4_WT - 0.39, y = 2.2, label = "sensitive"), color = "black", size = 5) +
  geom_segment(aes(x = mean_FoldChange_T4_WT - 3*sd_FoldChange_T3_WT, y = 2, xend = mean_FoldChange_T4_WT - sd_FoldChange_T4_WT, yend = 2), arrow = arrow(length = unit(0.4, "cm")), size = 1) +
  geom_segment(aes(x = mean_FoldChange_T4_WT - sd_FoldChange_T3_WT, y = 2, xend = mean_FoldChange_T4_WT - 3*sd_FoldChange_T4_WT, yend = 2), arrow = arrow(length = unit(0.4, "cm")), size = 1)
ggsave(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "/plots/T4_arrow_more.png"), width = 6, height = 4)


### Determine the thresholds based on graphs
## T3
limit_T3_resistant <- mean_FoldChange_T3_WT + 2*sd_FoldChange_T3_WT
limit_T3_hypersensitive <- mean_FoldChange_T3_WT - 2*sd_FoldChange_T3_WT

## T4
limit_T4_resistant <- mean_FoldChange_T4_WT + 2*sd_FoldChange_T4_WT
limit_T4_hypersensitive <- mean_FoldChange_T4_WT - 3*sd_FoldChange_T4_WT
limit_T4_moresensitive <- mean_FoldChange_T4_WT - sd_FoldChange_T4_WT

### Add information on resistance based on AUC, T3 and T4 foldchanges
## T3
c_temp_growth$T3_result <- with(c_temp_growth, 
                                    ifelse(FoldChange_gene_T3 > limit_T3_resistant, "resistant", 
                                           ifelse(FoldChange_gene_T3 < limit_T3_hypersensitive, "hypersensitive", "no effect")))
## T4
c_temp_growth$T4_result <- with(c_temp_growth, ifelse(FoldChange_gene_T4 > limit_T4_resistant, "resistant", 
                                                              ifelse(FoldChange_gene_T4 < limit_T4_hypersensitive, "hypersensitive", 
                                                                     ifelse(FoldChange_gene_T4 > limit_T4_hypersensitive & FoldChange_gene_T4 < limit_T4_moresensitive, "more sensitive", 
                                                                            "no effect"))))

## New limits for T3 resistant and more resistant
resistant_1sd <- mean_FoldChange_T3_WT + sd_FoldChange_T3_WT
resistant_3sd <- mean_FoldChange_T3_WT + 3*sd_FoldChange_T3_WT

# add a column with this new information
c_temp_growth$T3_more_less_resistant <- with(c_temp_growth, 
                                                 ifelse(FoldChange_gene_T3 > resistant_3sd, "3*sd resistant", 
                                                        ifelse(FoldChange_gene_T3 > resistant_1sd, "1*sd resistant", "no resistance")))


############################################################

########### Selection resistant strains
resistant_selection_T3 <- filter(c_temp_growth, T3_result == "resistant")
resistant_selection_T3$replicate <- duplicated(resistant_selection_T3$gene_name)
resistant_selection_T3 <- resistant_selection_T3 %>% drop_na(locus_tag)

########### Selection hypersensitive strains
hypersensitive_selection_T4 <- filter(c_temp_growth, T4_result == "hypersensitive")
hypersensitive_selection_T4$replicate <- duplicated(hypersensitive_selection_T4$gene_name)
hypersensitive_selection_T4 <- hypersensitive_selection_T4 %>% drop_na(locus_tag)

moresensitive_selection_T4 <- filter(c_temp_growth, T4_result == "more sensitive")
moresensitive_selection_T4$replicate <- duplicated(moresensitive_selection_T4$gene_name)
moresensitive_selection_T4 <- moresensitive_selection_T4 %>% drop_na(locus_tag)

############################################################

############################################################
### Add gene ontology (GO) information #####################
############################################################

#### Upload the data you need ####
# upload gpad file
url_gpad <- "http://current.geneontology.org/annotations/ecocyc.gpad.gz"
tmp_gpad <- tempfile()
download.file(url_gpad, tmp_gpad)
data_gpad <- read.csv(gzfile(tmp_gpad), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# upload gpi file
url_gpi <- "http://current.geneontology.org/annotations/ecocyc.gpi.gz"
tmp_gpi <- tempfile()
download.file(url_gpi, tmp_gpi)
data_gpi <- read.csv(gzfile(tmp_gpi), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# change column names
colnames_gpad <- c("DB", "DB_Object_ID", "Qualifier", "GO_ID", "DB:Reference", "Evidence_Code", "With_From", "Interacting_taxon_ID", "Date", "Assigned_by", "Annotation_extension", "Annotation_properties")
colnames(data_gpad) <- colnames_gpad
colnames_gpi <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "DB_Object_Name", "DB_Object_Synonym(s)", "DB_Object_Type", "Taxon", "Parent_Object_ID", "DB_Xref(s)", "Properties")
colnames(data_gpi) <- colnames_gpi
data_gpad <- data_gpad[-c(1), ]

# select only columns we need and rbind
gpad <- data_gpad[2:4]
gpi <- data_gpi[2:3]
go_terms <- right_join(gpad, gpi, by = "DB_Object_ID")
colnames(go_terms)[4] <- "gene_name"
go_terms_unique<- go_terms[!duplicated(go_terms),]

## create a dataset with the selected genes to be then able to make a loop
list_datasets <- list(resistant_selection_T3, hypersensitive_selection_T4, moresensitive_selection_T4)
names(list_datasets) <- c("resistant_selection_T3", "hypersensitive_selection_T4", "moresensitive_selection_T4")
datasets_withGO <- list()
datasets_withGO_count <- list()

### make a wide to long table
for(F in names(list_datasets)){
  print(F)
  test_table <- list_datasets[[F]]
  table_hits <- test_table[4:5]
  names(table_hits)[2] <- "gene_name"
  table_hits$gene_name <- as.character(table_hits$gene_name)
  table_hits <- table_hits[!duplicated(table_hits),]
  table_hits_go <- right_join(go_terms_unique, table_hits, by = "gene_name")
  table_hits_go$GO_name <- Term(table_hits_go$GO_ID)
  table_hits_go_count <- table_hits_go %>% 
    group_by(GO_ID, GO_name) %>% 
    tally() %>% 
    arrange(desc(n))
  table_hits_go_count$GO_ontology <- Ontology(table_hits_go_count$GO_ID)
  table_hits_go_count$GO_synonym <- Synonym(table_hits_go_count$GO_ID)
  table_hits_go_count[is.na(table_hits_go_count)] <- "NO ASSOCIATED GO TERM"
  table_hits_go_count <- table_hits_go_count[c("GO_ID", "GO_name", "GO_ontology", "GO_synonym", "n")]
  datasets_withGO_count[[F]] <- table_hits_go_count
  datasets_withGO[[F]] <- table_hits_go
}

### extract the different tables from the dataset
## data with GO count
resistant_selection_T3_GO_count <- datasets_withGO_count[["resistant_selection_T3"]]
hypersensitive_selection_T4_GO_count <- datasets_withGO_count[["hypersensitive_selection_T4"]]
moresensitive_selection_T4_GO_count <- datasets_withGO_count[["moresensitive_selection_T4"]]

## data with GO ID and GO name
resistant_selection_T3_GO <- datasets_withGO[["resistant_selection_T3"]]
hypersensitive_selection_T4_GO <- datasets_withGO[["hypersensitive_selection_T4"]]
moresensitive_selection_T4_GO <- datasets_withGO[["moresensitive_selection_T4"]]

############################################################
### GO term analysis for hypersensitive strains ############

#add GO ontology
hypersensitive_selection_T4_GO$GO_ontology <- Ontology(hypersensitive_selection_T4_GO$GO_ID)
hypersensitive_selection_T4_GO <- hypersensitive_selection_T4_GO[!duplicated(hypersensitive_selection_T4_GO),]

#make 3 dataframes based on ontology
hypersensitive_selection_T4_GO_BP <- filter(hypersensitive_selection_T4_GO, GO_ontology == "BP")
hypersensitive_selection_T4_GO_CC <- filter(hypersensitive_selection_T4_GO, GO_ontology == "CC")
hypersensitive_selection_T4_GO_MF <- filter(hypersensitive_selection_T4_GO, GO_ontology == "MF")

#count for each the number of genes
hypersensitive_selection_T4_GO_BP_count <- hypersensitive_selection_T4_GO_BP %>% 
  group_by(GO_ID, GO_name) %>%
  tally() %>%
  arrange(desc(n))
hypersensitive_selection_T4_GO_CC_count <- hypersensitive_selection_T4_GO_CC %>% 
  group_by(GO_ID, GO_name) %>%
  tally() %>%
  arrange(desc(n))
hypersensitive_selection_T4_GO_MF_count <- hypersensitive_selection_T4_GO_MF %>% 
  group_by(GO_ID, GO_name) %>%
  tally() %>%
  arrange(desc(n))

# make a plot for BP 
hypersensitive_selection_T4_GO_BP_count_plot <- ungroup(hypersensitive_selection_T4_GO_BP_count) %>% top_n(n=10, n)
g <- ggplot(hypersensitive_selection_T4_GO_BP_count_plot, aes(x = reorder(GO_name, n), y = n)) 
g + geom_bar(stat = "identity", fill = "plum1", color = "black", width = 0.7, size = 0.3) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw() +
  labs(x = NULL, y = "Number of hypersensitive strains", title = "Biological processes (BP)")
ggsave(paste0("./BP_plot.png"), width = 4, height = 4)

# make a plot for CC 
hypersensitive_selection_T4_GO_CC_count_plot <- ungroup(hypersensitive_selection_T4_GO_CC_count) %>% top_n(n=10, n)
g <- ggplot(hypersensitive_selection_T4_GO_CC_count_plot, aes(x = reorder(GO_name, n), y = n)) 
g + geom_bar(stat = "identity", fill = "plum2", color = "black", width = 0.7, size = 0.3) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw() +
  labs(x = NULL, y = "Number of hypersensitive strains", title = "Cellular components (CC)")
ggsave(paste0("./CC_plot.png"), width = 4, height = 4)

# make a plot for MF 
hypersensitive_selection_T4_GO_MF_count_plot <- ungroup(hypersensitive_selection_T4_GO_MF_count) %>% top_n(n=10, n)
g <- ggplot(hypersensitive_selection_T4_GO_MF_count_plot, aes(x = reorder(GO_name, n), y = n)) 
g + geom_bar(stat = "identity", fill = "plum3", color = "black", width = 0.7, size = 0.3) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw() +
  labs(x = NULL, y = "Number of hypersensitive strains", title = "Molecular function (MF)")
ggsave(paste0("./MF_plot.png"), width = 4, height = 4)

############################################################
############ Add KEGG information ##########################
############################################################

##libraries
library(KEGGREST)
library(reshape2)

## upload KEIO library list
## select only locus_tag and add "eco:" to have the Kegg ID
KEIO_samples <- read_excel("c://Users/theinone/Desktop/Screening_data/KEIO_samples.xlsx")
KEIO_samples_analysis <- KEIO_samples %>% 
  dplyr::select(locus_tag, gene_name) %>% 
  drop_na() %>% 
  unique()
KEIO_samples$eco_ID <- paste("eco", KEIO_samples$locus_tag, sep = ":")

## assign path_ID and description
# extract information in the form of a "named vector"
ko <- keggLink("ko", "eco")
pathway_ID <- keggLink("pathway", "eco")
pathway_description <- keggList("pathway", "eco")
# convert to dataframe
df_ko <- data.frame(eco_ID=names(ko), ko_ID=ko, row.names = NULL, stringsAsFactors = FALSE)
df_pathways <- data.frame(eco_ID=names(pathway_ID), pathway_ID=pathway_ID, row.names = NULL, stringsAsFactors = FALSE)
df_description <- data.frame(pathway_ID=names(pathway_description), pathway_description=pathway_description, row.names = NULL, stringsAsFactors = FALSE)
#merge tables
KEGG_df_temp <- left_join(df_ko, df_pathways, by = "eco_ID")
KEGG_df <- left_join(KEGG_df_temp, df_description, by ="pathway_ID")
KEIO_samples_ko <- left_join(KEIO_samples, KEGG_df, by = "eco_ID")

##################################################
########## COUNT NUMBER GENES / PATHWAY ##########
##################################################

pathways_KEIO <- KEIO_samples_ko %>% 
  dplyr::select(c("locus_tag", "gene_name", "pathway_ID", "pathway_description"))
pathways_KEIO <- pathways_KEIO[!is.na(pathways_KEIO$locus_tag),]
pathways_KEIO_unique <- pathways_KEIO[!duplicated(pathways_KEIO),]
pathways_KEIO_count <- pathways_KEIO_unique %>% 
  group_by(pathway_ID, pathway_description) %>%
  tally() %>%
  arrange(desc(n))
# rename column to add that it's in the whole genome
names(pathways_KEIO_count)[names(pathways_KEIO_count) == "n"] <- "number_genome"
pathways_KEIO_count[is.na(pathways_KEIO_count)] <- "NO ASSOCIATED KEGG PATHWAY"

## number of unique genes in the KEIO collection
unique_genes <- KEIO_samples %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_unique_genes <- as.numeric(nrow(unique_genes))

## number of unique genes in the KEIO collection genes associated to a KEGG pathway
unique_genes_pathway <- pathways_KEIO_unique[!is.na(pathways_KEIO_unique$pathway_ID),]
unique_genes_pathway <- unique_genes_pathway %>% 
  dplyr::select(locus_tag) %>% 
  unique()
number_unique_genes_pathway <- as.numeric(nrow(unique_genes_pathway))

## calculate the fraction for each pathway
pathways_KEIO_count$fraction_genome <- NA
pathways_KEIO_count$fraction_genome <- round((pathways_KEIO_count$number_genome / number_unique_genes), digits = 3)

##########################################################
########## MERGE KEGG INFO WITH SELECTION GENES ##########
##########################################################

## open the list with all information about the genes and pathways 
KEIO_KEGG_list <- KEIO_samples_ko[c("locus_tag", "eco_ID", "pathway_ID", "pathway_description")]
KEIO_KEGG_percentage_genome <- pathways_KEIO_count %>% dplyr::select(pathway_ID, number_genome, fraction_genome)

#### Resistant T3 selection #####
## add KEGG information to resistant strains (T3)
resistant_selection_T3_KEGG <- left_join(resistant_selection_T3[2:3], KEIO_KEGG_list, by = "locus_tag")
resistant_selection_T3_KEGG <- unique(resistant_selection_T3_KEGG)
count_resistant_selection_T3_KEGG <- resistant_selection_T3_KEGG %>%
  group_by(pathway_ID, pathway_description) %>%
  tally() %>%
  arrange(desc(n))
names(count_resistant_selection_T3_KEGG)[names(count_resistant_selection_T3_KEGG) == "n"] <- "number_T3"
count_resistant_selection_T3_KEGG[is.na(count_resistant_selection_T3_KEGG)] <- "NO ASSOCIATED KEGG PATHWAY"

## number of resistant genes (T3)
resistant_T3 <- resistant_selection_T3 %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_resistant_T3 <- as.numeric(nrow(resistant_T3))

## number of resistant genes associated to pathways (T3)
resistant_T3_pathways <- resistant_selection_T3_KEGG[!is.na(resistant_selection_T3_KEGG$pathway_ID),]
resistant_T3_pathways <- resistant_T3_pathways %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_resistant_T3_pathways <- as.numeric(nrow(resistant_T3_pathways))

## calculate the fraction for each pathway
count_resistant_selection_T3_KEGG$fraction_inselection <- NA
count_resistant_selection_T3_KEGG$fraction_inselection <- round((count_resistant_selection_T3_KEGG$number_T3 / number_resistant_T3), digits = 3)
count_resistant_selection_T3_KEGG <- left_join(count_resistant_selection_T3_KEGG, KEIO_KEGG_percentage_genome, by = "pathway_ID")

## calculate the fraction compared to the whole genome
count_resistant_selection_T3_KEGG$fraction_ingenome <- NA
count_resistant_selection_T3_KEGG$fraction_ingenome <- with(count_resistant_selection_T3_KEGG, round((number_T3 / number_genome), digits = 3))

## save as excel
write_xlsx(count_resistant_selection_T3_KEGG, file.path(paste0("./Results/", date, "/resistant/resistant_selection_T3_KEGG_count.xlsx")))
write_xlsx(resistant_selection_T3_KEGG, file.path(paste0("./Results/", date, "/resistant/resistant_selection_T3_KEGG_terms.xlsx")))     
#################################

#### Hypersensitive T4 selection #####
## add KEGG information to hypersensitive genes (T4)
hypersensitive_selection_T4_KEGG <- left_join(hypersensitive_selection_T4[2:3], KEIO_KEGG_list, by = "locus_tag")
hypersensitive_selection_T4_KEGG <- unique(hypersensitive_selection_T4_KEGG)
count_hypersensitive_selection_T4_KEGG <- hypersensitive_selection_T4_KEGG %>%
  group_by(pathway_ID, pathway_description) %>%
  tally() %>%
  arrange(desc(n))
names(count_hypersensitive_selection_T4_KEGG)[names(count_hypersensitive_selection_T4_KEGG) == "n"] <- "number_T4"
count_hypersensitive_selection_T4_KEGG[is.na(count_hypersensitive_selection_T4_KEGG)] <- "NO ASSOCIATED KEGG PATHWAY"

## number of hypersensitive genes (T4)
hypersensitive_T4 <- hypersensitive_selection_T4 %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_hypersensitive_T4 <- as.numeric(nrow(hypersensitive_T4))

## number of hypersensitive genes assigned to pathways(T4)
hypersensitive_genes_pathways <- hypersensitive_selection_T4_KEGG[!is.na(hypersensitive_selection_T4_KEGG$pathway_ID),]
hypersensitive_genes_pathways <- hypersensitive_genes_pathways %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_hypersensitive_T4_pathways <- as.numeric(nrow(hypersensitive_genes_pathways))

## calculate the fraction for each pathway
count_hypersensitive_selection_T4_KEGG$fraction_inselection <- NA
count_hypersensitive_selection_T4_KEGG$fraction_inselection <- round((count_hypersensitive_selection_T4_KEGG$number_T4 / number_hypersensitive_T4), digits = 3)
count_hypersensitive_selection_T4_KEGG <- left_join(count_hypersensitive_selection_T4_KEGG, KEIO_KEGG_percentage_genome, by = "pathway_ID")

## calculate the fraction compared to the whole genome
count_hypersensitive_selection_T4_KEGG$fraction_ingenome <- NA
count_hypersensitive_selection_T4_KEGG$fraction_ingenome <- with(count_hypersensitive_selection_T4_KEGG, round((number_T4 / number_genome), digits = 3))

## save as excel
write_xlsx(count_hypersensitive_selection_T4_KEGG, file.path(paste0("./Results/", date, "/hypersensitive/hypersensitive_selection_T4_KEGG_count.xlsx")))
write_xlsx(hypersensitive_selection_T4_KEGG, file.path(paste0("./Results/", date, "/hypersensitive/hypersensitive_selection_T4_KEGG_terms.xlsx")))
#################################


#### Moresensitive T4 selection #####
## add KEGG information to moresensitive genes (T4)
moresensitive_selection_T4_KEGG <- left_join(moresensitive_selection_T4[2:3], KEIO_KEGG_list, by = "locus_tag")
moresensitive_selection_T4_KEGG <- unique(moresensitive_selection_T4_KEGG)
count_moresensitive_selection_T4_KEGG <- moresensitive_selection_T4_KEGG %>%
  group_by(pathway_ID, pathway_description) %>%
  tally() %>%
  arrange(desc(n))
names(count_moresensitive_selection_T4_KEGG)[names(count_moresensitive_selection_T4_KEGG) == "n"] <- "number_T4"
count_moresensitive_selection_T4_KEGG[is.na(count_moresensitive_selection_T4_KEGG)] <- "NO ASSOCIATED KEGG PATHWAY"

## number of moresensitive genes (T4)
moresensitive_T4 <- moresensitive_selection_T4 %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_moresensitive_T4 <- as.numeric(nrow(moresensitive_T4))

## number of moresensitive genes assigned to pathways(T4)
moresensitive_genes_pathways <- moresensitive_selection_T4_KEGG[!is.na(moresensitive_selection_T4_KEGG$pathway_ID),]
moresensitive_genes_pathways <- moresensitive_genes_pathways %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_moresensitive_T4_pathways <- as.numeric(nrow(moresensitive_genes_pathways))

## calculate the percentage for each pathway
count_moresensitive_selection_T4_KEGG$fraction_inselection <- NA
count_moresensitive_selection_T4_KEGG$fraction_inselection <- round((count_moresensitive_selection_T4_KEGG$number_T4 / number_moresensitive_T4), digits = 2)
count_moresensitive_selection_T4_KEGG <- left_join(count_moresensitive_selection_T4_KEGG, KEIO_KEGG_percentage_genome, by = "pathway_ID")

## calculate the percentage compared to the whole genome
count_moresensitive_selection_T4_KEGG$fraction_ingenome <- NA
count_moresensitive_selection_T4_KEGG$fraction_ingenome <- with(count_moresensitive_selection_T4_KEGG, round((number_T4 / number_genome), digits = 2))

## save as excel
write_xlsx(count_moresensitive_selection_T4_KEGG, file.path(paste0("./Results/", date, "/hypersensitive/moresensitive_selection_T4_KEGG_count.xlsx")))
write_xlsx(moresensitive_selection_T4_KEGG, file.path(paste0("./Results/", date, "/hypersensitive/moresensitive_selection_T4_KEGG_terms.xlsx")))
#################################

## save information on number of genes in each set
number_genes <- data.frame(number_resistant_T3, number_resistant_T3_pathways, number_hypersensitive_T4, number_hypersensitive_T4_pathways, number_moresensitive_T4, number_moresensitive_T4_pathways)
number_genes_transposed <- as.data.frame(t(number_genes))
number_genes_transposed$selection <- c("resistant T3", "resistant T3 in pathways", "hypersensitive T4", "hypersensitive T4 in pathways", "moresensitive T4", "moresensitive T4 in pathways")
names(number_genes_transposed)[1] <- "number"
number_genes_transposed <- number_genes_transposed[2:1]
# add the number of genes in the KEIO collection
number_genes_transposed <- number_genes_transposed %>% add_row(selection = "KEIO", number = 3896) %>% add_row(selection = "KEIO in pathways", number = 1362)
write_xlsx(x = number_genes_transposed, file.path(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "/number_genes.xlsx")))


########################################################
############## plots with KEGG information #############
########################################################

## upload genes sets 
## /!\ check that you have the latest version /!\ ##
resistant_T3 <- count_resistant_selection_T3_KEGG
hypersensitive <- count_hypersensitive_selection_T4_KEGG
#hypersensitive <- read_excel("c://Users/theinone/Desktop/Screening_data/Results/2020_07_16/hypersensitive/hypersensitive_selection_T4_KEGG_count.xlsx")
#resistant_AUC <- read_excel("c://Users/theinone/Desktop/Screening_data/Results/2020_07_16/resistant/resistant_selection_AUC_KEGG_count.xlsx")

## make a dataset with the two gene sets
datasets <- list(hypersensitive, resistant_T3)
names(datasets) <- c("hypersensitive", "resistant_T3")
datasets_selection <- list()

### keep only pathways that have at least 5 genes and add information about each gene set (hypersensitive or resistant)
for(F in names(datasets)){
  print(F)
  input_table <- datasets[[F]]
  selection <- filter(input_table, number_genome > 4)
  selection$pathway_description <- sub(pattern = "- Escherichia coli K-12 MG1655", "", selection$pathway_description) 
  selection <- selection %>% dplyr::select(pathway_ID, pathway_description, fraction_inselection, fraction_ingenome, contains("number"))
  colnames(selection)[3] <- paste(F, colnames(selection[,c(3)]), sep="_")
  colnames(selection)[4] <- paste(F, "fraction_ingenome", sep="_")
  datasets_selection[[F]] <- selection
}

## retrieve dataframes
hypersensitive_selection <- datasets_selection[["hypersensitive"]]
resistant_T3_selection <- datasets_selection[["resistant_T3"]]

## make one dataframe with all information for hypersensitive genes
all_selection <- full_join(hypersensitive_selection, resistant_T3_selection, by = "pathway_description")
all_selection <- ungroup(all_selection)

## make three dataframes
# one with information about percentage in selection
selection_in_selection <- all_selection %>% dplyr::select(pathway_description, contains("fraction_inselection"))
long_selection_in_selection <- melt(selection_in_selection, 
                                    variable.name = "phenotype", 
                                    value.name = "fraction_in_selection", 
                                    id.vars = c("pathway_description"))
long_selection_in_selection$phenotype <- sub(pattern = "_fraction_inselection", "", long_selection_in_selection$phenotype)
long_selection_in_selection$fraction_in_selection[is.na(long_selection_in_selection$fraction_in_selection)] <- 0
long_selection_in_selection <- long_selection_in_selection[!grepl("NO ASSOCIATED", long_selection_in_selection$pathway_description),]

# one with information about percentage in genome
selection_in_genome <- all_selection %>% dplyr::select(pathway_description, contains("fraction_ingenome"))
long_selection_in_genome <- melt(selection_in_genome, 
                                 variable.name = "phenotype", 
                                 value.name = "fraction_in_genome", 
                                 id.vars = c("pathway_description"))
long_selection_in_genome$phenotype <- sub(pattern = "_fraction_ingenome", "", long_selection_in_genome$phenotype)
long_selection_in_genome$fraction_in_genome[is.na(long_selection_in_genome$fraction_in_genome)] <- 0
long_selection_in_genome <- long_selection_in_genome[!grepl("NO ASSOCIATED", long_selection_in_genome$pathway_description),]



##################################################
################## MAKE GRAPHS ###################
##################################################

library(dplyr)
library(ggplot2)

## point graphs with all data
color_selection <- c("hypersensitive" = "lightskyblue", "resistant_AUC" = "pink", "resistant_T3" = "plum")
g <- ggplot(long_selection_in_genome, aes(x = reorder(pathway_description, fraction_in_genome), y = fraction_in_genome))
g + geom_point(aes(color = phenotype)) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Fraction" , title = "Fraction in all KEIO genes") +
  scale_color_manual(values = color_selection)

p <- ggplot(long_selection_in_selection, aes(x = reorder(pathway_description, fraction_in_selection), y = fraction_in_selection))
p + geom_point(aes(color = phenotype)) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Fraction", title = "Fraction in selected genes (omit genes not assigned to pathways)") +
  scale_color_manual(values = color_selection )

## point graphs with only hypersensitive genes
hypersensitive_selection_plot <- ungroup(hypersensitive_selection) %>% filter(number_T4 > 1) %>%
  top_n(n=10, hypersensitive_fraction_ingenome)

ggplot(hypersensitive_selection_plot, aes(x=reorder(pathway_description, hypersensitive_fraction_ingenome), y = hypersensitive_fraction_ingenome)) +
  geom_point(shape = 21, color = "black", fill = "pink", aes(size=number_T4)) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Fraction among KEIO library strains", title = "Hypersensitive deletion strains", size = "Number strains") +
  scale_y_continuous(limits = c(0.2, 0.6))
ggsave(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "_hypersensitive_genes.png"), width = 5, height = 4)

## point graphs with only resisant genes
resistant_T3_selection_plot <- ungroup(resistant_T3_selection) %>% 
  filter(number_T3 > 1) %>%
  top_n(n=10, resistant_T3_fraction_ingenome)

ggplot(resistant_T3_selection_plot, aes(x=reorder(pathway_description, resistant_T3_fraction_ingenome), y = resistant_T3_fraction_ingenome)) +
  geom_point(shape = 21, color = "black", fill = "blue", aes(size=number_T3)) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Fraction among KEIO library strains", title = "Resistant deletion strains T3", size = "Number strains") 
ggsave(paste0("c://Users/theinone/Desktop/Screening_data/Results/", date, "_resistant_AUC_genes.png"), width = 5, height = 4)

##################################################
################### END script ###################
##################################################

