################################################################################
############################### GO term analysis ###############################
################################################################################

## libraries
library(tidyverse)
library(readxl)
library(tools)
library(stringr)
library(GO.db)

### upload the data you need 
# upload gpad file from GO database
url_gpad <- "http://current.geneontology.org/annotations/ecocyc.gpad.gz"
tmp_gpad <- tempfile()
download.file(url_gpad, tmp_gpad)
data_gpad <- read.csv(gzfile(tmp_gpad), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# upload gpi file from GIO database
url_gpi <- "http://current.geneontology.org/annotations/ecocyc.gpi.gz"
tmp_gpi <- tempfile()
download.file(url_gpi, tmp_gpi)
data_gpi <- read.csv(gzfile(tmp_gpi), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# change column names for each file
colnames_gpad <- c("DB", "DB_Object_ID", "Qualifier", "GO_ID", "DB:Reference", "Evidence_Code", "With_From", "Interacting_taxon_ID", "Date", "Assigned_by", "Annotation_extension", "Annotation_properties")
colnames(data_gpad) <- colnames_gpad
colnames_gpi <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "DB_Object_Name", "DB_Object_Synonym(s)", "DB_Object_Type", "Taxon", "Parent_Object_ID", "DB_Xref(s)", "Properties")
colnames(data_gpi) <- colnames_gpi
data_gpad <- data_gpad[-c(1), ]

# select only columns you need and rbind
gpad <- data_gpad[2:4]
gpi <- data_gpi[2:3]
go_terms <- right_join(gpad, gpi, by = "DB_Object_ID")
colnames(go_terms)[4] <- "gene_name"
go_terms_unique<- go_terms[!duplicated(go_terms),]

### add the GO information to the set of hypersensitive genes

# select only the column with locus_tag and gene name
table_hits <- hypersensitive_selection_T4[2:3]
table_hits$gene_name <- as.character(table_hits$gene_name)
table_hits <- table_hits[!duplicated(table_hits),]

# add the GO annotation
hypersensitive_selection_T4_GO <- right_join(go_terms_unique, table_hits, by = "gene_name")

# obtain the GO name
hypersensitive_selection_T4_GO$GO_name <- Term(hypersensitive_selection_T4_GO$GO_ID)

# add the GO ontolology information
hypersensitive_selection_T4_GO$GO_ontology <- Ontology(hypersensitive_selection_T4_GO$GO_ID)
hypersensitive_selection_T4_GO <- hypersensitive_selection_T4_GO[!duplicated(hypersensitive_selection_T4_GO),]

# divide in 3 datasets, one for each ontology
hypersensitive_selection_T4_GO_BP <- filter(hypersensitive_selection_T4_GO, GO_ontology == "BP")
hypersensitive_selection_T4_GO_CC <- filter(hypersensitive_selection_T4_GO, GO_ontology == "CC")
hypersensitive_selection_T4_GO_MF <- filter(hypersensitive_selection_T4_GO, GO_ontology == "MF")

# count the number of genes for each GO term in each GO ontology
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
  
## make a plot for the top 10 GO terms in each GO ontology
# biological processes (BP)
hypersensitive_selection_T4_GO_BP_count_plot <- ungroup(hypersensitive_selection_T4_GO_BP_count) %>% top_n(n=10, n)
g <- ggplot(hypersensitive_selection_T4_GO_BP_count_plot, aes(x = reorder(GO_name, n), y = n)) 
g + geom_bar(stat = "identity", fill = "plum1", color = "black", width = 0.7, size = 0.3) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw() +
  labs(x = NULL, y = "Number of hypersensitive strains", title = "Biological processes (BP)")

# cellular component (CC)
hypersensitive_selection_T4_GO_CC_count_plot <- ungroup(hypersensitive_selection_T4_GO_CC_count) %>% top_n(n=10, n)
g <- ggplot(hypersensitive_selection_T4_GO_CC_count_plot, aes(x = reorder(GO_name, n), y = n)) 
g + geom_bar(stat = "identity", fill = "plum2", color = "black", width = 0.7, size = 0.3) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw() +
  labs(x = NULL, y = "Number of hypersensitive strains", title = "Cellular components (CC)")

# molecular function (MF)
hypersensitive_selection_T4_GO_MF_count_plot <- ungroup(hypersensitive_selection_T4_GO_MF_count) %>% top_n(n=10, n)
g <- ggplot(hypersensitive_selection_T4_GO_MF_count_plot, aes(x = reorder(GO_name, n), y = n)) 
g + geom_bar(stat = "identity", fill = "plum3", color = "black", width = 0.7, size = 0.3) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw() +
  labs(x = NULL, y = "Number of hypersensitive strains", title = "Molecular function (MF)")

#################################################################################
########################## END script GO term analysis ##########################
#################################################################################
