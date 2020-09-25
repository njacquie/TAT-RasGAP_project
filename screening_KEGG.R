##libraries
library(KEGGREST)
library(tidyverse)
library(reshape2)

### add KEGG information to the E. coli KEIO collection

# upload the E. coli Keio collection - Plate maps (available online)
KEIO_samples <- read_excel("./KEIO_samples.xlsx")

# select only locus_tag (also called E. coli MG1655 B id) and add "eco:" to have the Kegg ID
KEIO_samples_analysis <- KEIO_samples %>% 
  dplyr::select(locus_tag, gene_name) %>% 
  drop_na() %>% 
  unique()
KEIO_samples$eco_ID <- paste("eco", KEIO_samples$locus_tag, sep = ":")

# extract path_ID and description information in the form of a "named vector"
ko <- keggLink("ko", "eco")
pathway_ID <- keggLink("pathway", "eco")
pathway_description <- keggList("pathway", "eco")

# convert to dataframe
df_ko <- data.frame(eco_ID=names(ko), ko_ID=ko, row.names = NULL, stringsAsFactors = FALSE)
df_pathways <- data.frame(eco_ID=names(pathway_ID), pathway_ID=pathway_ID, row.names = NULL, stringsAsFactors = FALSE)
df_description <- data.frame(pathway_ID=names(pathway_description), pathway_description=pathway_description, row.names = NULL, stringsAsFactors = FALSE)

# merge tables
KEGG_df_temp <- left_join(df_ko, df_pathways, by = "eco_ID")
KEGG_df <- left_join(KEGG_df_temp, df_description, by ="pathway_ID")
KEIO_samples_ko <- left_join(KEIO_samples, KEGG_df, by = "eco_ID")

### count the number of gene of the E. coli KEIO collection per pathway

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

# number of unique genes in the KEIO collection
unique_genes <- KEIO_samples %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_unique_genes <- as.numeric(nrow(unique_genes))

# number of unique genes in the KEIO collection genes associated to a KEGG pathway
unique_genes_pathway <- pathways_KEIO_unique[!is.na(pathways_KEIO_unique$pathway_ID),]
unique_genes_pathway <- unique_genes_pathway %>% 
  dplyr::select(locus_tag) %>% 
  unique()
number_unique_genes_pathway <- as.numeric(nrow(unique_genes_pathway))

# calculate the fraction for each pathway
pathways_KEIO_count$fraction_genome <- NA
pathways_KEIO_count$fraction_genome <- round((pathways_KEIO_count$number_genome / number_unique_genes), digits = 3)

### count number of gene deletion inducing hypersensitivity in each KEGG pathway

hypersensitive_selection_T4_KEGG <- left_join(hypersensitive_selection_T4[2:3], KEIO_KEGG_list, by = "locus_tag")
hypersensitive_selection_T4_KEGG <- unique(hypersensitive_selection_T4_KEGG)
count_hypersensitive_selection_T4_KEGG <- hypersensitive_selection_T4_KEGG %>%
  group_by(pathway_ID, pathway_description) %>%
  tally() %>%
  arrange(desc(n))
names(count_hypersensitive_selection_T4_KEGG)[names(count_hypersensitive_selection_T4_KEGG) == "n"] <- "number_T4"
count_hypersensitive_selection_T4_KEGG[is.na(count_hypersensitive_selection_T4_KEGG)] <- "NO ASSOCIATED KEGG PATHWAY"

# number of hypersensitive genes (T4)
hypersensitive_T4 <- hypersensitive_selection_T4 %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_hypersensitive_T4 <- as.numeric(nrow(hypersensitive_T4))

# number of hypersensitive genes assigned to pathways(T4)
hypersensitive_genes_pathways <- hypersensitive_selection_T4_KEGG[!is.na(hypersensitive_selection_T4_KEGG$pathway_ID),]
hypersensitive_genes_pathways <- hypersensitive_genes_pathways %>% 
  dplyr::select(locus_tag) %>% 
  drop_na() %>%
  unique()
number_hypersensitive_T4_pathways <- as.numeric(nrow(hypersensitive_genes_pathways))

# calculate the fraction for each pathway
count_hypersensitive_selection_T4_KEGG$fraction_inselection <- NA
count_hypersensitive_selection_T4_KEGG$fraction_inselection <- round((count_hypersensitive_selection_T4_KEGG$number_T4 / number_hypersensitive_T4), digits = 3)
count_hypersensitive_selection_T4_KEGG <- left_join(count_hypersensitive_selection_T4_KEGG, KEIO_KEGG_percentage_genome, by = "pathway_ID")

# calculate the fraction compared to the whole genome
count_hypersensitive_selection_T4_KEGG$fraction_ingenome <- NA
count_hypersensitive_selection_T4_KEGG$fraction_ingenome <- with(count_hypersensitive_selection_T4_KEGG, round((number_T4 / number_genome), digits = 3))

### make plot

# select pathways that have at least 5 genes and add information 
hypersensitive_selection <- filter(count_hypersensitive_selection_T4_KEGG, number_genome > 4)
hypersensitive_selection$pathway_description <- sub(pattern = "- Escherichia coli K-12 MG1655", "", hypersensitive_selection$pathway_description)
hypersensitive_selection <- hypersensitive_selection %>% 
  dplyr::select(pathway_ID, pathway_description, fraction_inselection, fraction_ingenome, contains("number"))

# select 
hypersensitive_selection_plot <- ungroup(hypersensitive_selection) %>% filter(number_T4 > 1) %>%
  top_n(n=10, hypersensitive_fraction_ingenome)

ggplot(hypersensitive_selection_plot, aes(x=reorder(pathway_description, hypersensitive_fraction_ingenome), y = hypersensitive_fraction_ingenome)) +
  geom_point(shape = 21, color = "black", fill = "pink", aes(size=number_T4)) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Fraction among KEIO library strains", title = "Hypersensitive deletion strains", size = "Number strains") +
  scale_y_continuous(limits = c(0.2, 0.6))




