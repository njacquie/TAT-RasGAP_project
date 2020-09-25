## libraries
library(tidyverse)
library(readxl)
library(writexl)
library(tools)
library(stringr)

## upload dataset
c_temp <- read_excel("./dataset1.xlsx")

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

## calculate mean and sd at T3 and T4
## T3
mean_FoldChange_T3_WT <- mean(c_temp$FoldChange_WT_T3)
sd_FoldChange_T3_WT <- sd(c_temp$FoldChange_WT_T3)
maxValue_T3 <- max(c_temp$FoldChange_gene_T3)

## T4
mean_FoldChange_T4_WT <- mean(new_c_temp$FoldChange_WT_T4)
sd_FoldChange_T4_WT <- sd(new_c_temp$FoldChange_WT_T4)

## select only genes that did grow more to more than 0.4 OD600 at 24h without peptide
new_c_temp_growth <- filter(c_temp, gene_CTR_T4_minusT0 > 0.4)

### plot the distribution of the points to select the best thresholds
## T3 density plot with arrow
p <- ggplot(new_c_temp_growth, aes(FoldChange_gene_T3))
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

## T4 density plot with arrows
q <- ggplot(new_c_temp_growth, aes(FoldChange_gene_T4))
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

### Add information on resistance based on T3 and T4 foldchanges
## T3
new_c_temp_growth$T3_result <- with(new_c_temp_growth, 
                                    ifelse(FoldChange_gene_T3 > limit_T3_resistant, "resistant", 
                                           ifelse(FoldChange_gene_T3 < limit_T3_hypersensitive, "hypersensitive", "no effect")))
## T4
new_c_temp_growth$T4_result <- with(new_c_temp_growth, ifelse(FoldChange_gene_T4 > limit_T4_resistant, "resistant", 
                                                              ifelse(FoldChange_gene_T4 < limit_T4_hypersensitive, "hypersensitive", 
                                                                     ifelse(FoldChange_gene_T4 > limit_T4_hypersensitive & FoldChange_gene_T4 < limit_T4_moresensitive, "more sensitive", 
                                                                            "no effect"))))

### Make datasets with resistant and hypersensitive strains
# resistant strains
resistant_selection_T3 <- filter(new_c_temp_growth, T3_result == "resistant")
resistant_selection_T3$replicate <- duplicated(resistant_selection_T3$gene_name)
resistant_selection_T3 <- resistant_selection_T3 %>% drop_na(locus_tag)

# hypersensitive strains
hypersensitive_selection_T4 <- filter(new_c_temp_growth, T4_result == "hypersensitive")
hypersensitive_selection_T4$replicate <- duplicated(hypersensitive_selection_T4$gene_name)
hypersensitive_selection_T4 <- hypersensitive_selection_T4 %>% drop_na(locus_tag)

#################################################################################################
######## From here, use hypersensitive dataset for GO and KEGG analyses #########################
