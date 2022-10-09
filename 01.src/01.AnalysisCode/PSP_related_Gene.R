###########################################################
#      Prj: Multi-omics data analysis of P.sibiricum.
#      Assignment: PSP genes
#      Date: Mar 21, 2022
#      Author: Shawn Wang <shawnwang2016@126.com>
#      Location: HENU, Kaifeng, Henan, China
###########################################################
dir.create("../03.Progress/03.Transcriptome/PSP")
if (getwd()!="/Users/shawn/My_Repo/P.sibiricum/01.src") {
  setwd("~/My_Repo/P.sibiricum/03.Progress/03.Transcriptome/PSP/")
} else {
  setwd("../03.Progress/03.Transcriptome/PSP")
}


# packages ----------------------------------------------------------------

library(tidyverse)


psp_gene_anno <- read.delim("../PSP_related_Gene.txt",header = T)
# expMat_clean_merge <-
#   expMat_clean %>%
#   pivot_longer(
#     !gene_id,values_to = 'value',names_to = 'sample'
#   ) %>%
#   mutate(sample = gsub("..$","",sample)) %>%
#   group_by(gene_id,sample) %>%
#   summarise(fpkm = sum(value)) %>%
#   pivot_wider(
#     names_from = sample,values_from = fpkm
#   )
# writexl::write_xlsx(
#  x =  list(
#     split = expMat_clean,
#     merge = expMat_clean_merge
#   ),
#  path = "../FPKM.xlsx"
# )
DETs <- readxl::read_xlsx("../../02.Pacbio/Exp_data_split/DET/FPKM_merge_all.xlsx")

# merge -------------------------------------------------------------------

psp_mat <-
  inner_join(psp_gene_anno,DETs,c("GeneID" = "Gene_ID"))

psp_mat_clean <-
  psp_mat %>%
  group_by(category) %>%
  slice_min(order_by = E_value,n = 5) %>%
  ungroup() %>%
  select(GeneID,category,contains('Size')) %>%
  filter(rowSums(.[-(1:2)] > 1) >2) %>%
  mutate(
    New_ID = paste0(
      "HJ_iso_",
      str_pad(
        string = gsub("HJ_1-10k_transcript/","",GeneID),width = 6,side = 'left',pad = '0'
        )
    )
  ) %>% select(-GeneID)


library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
ht_mat <-
  psp_mat_clean %>%
  column_to_rownames("New_ID") %>%
  select(-category) %>%
  t() %>%
  scale() %>%
  t()

l_anno = rowAnnotation(
  tag = data.frame(
    row.names = psp_mat_clean$New_ID,
    category = psp_mat_clean$category
  ) %>% as.matrix()
)

ht_mat_plt <-
Heatmap(ht_mat,show_row_names = F,cluster_columns = F,left_annotation = l_anno,row_km = 6,row_km_repeats = 3000)

ht_mat_plt2 <- draw(ht_mat_plt)
ht_shiny(ht_mat_plt)




c6_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`6`,]

c4_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`4`,]

c3_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`3`,]



# sum ---------------------------------------------------------------------


psp_mat_sum =
  psp_mat_clean %>%
  select(-New_ID) %>%
  pivot_longer(
    !category,names_to = 'sample',values_to = 'value'
  ) %>%
  mutate(
    sample = gsub("..$","",sample)
  ) %>%
  group_by(category,sample) %>%
  summarise(sum = sum(value)) %>%
  pivot_wider(
    names_from = sample,values_from = sum
  )
ht_sum <-
  psp_mat_sum %>%
  column_to_rownames("category") %>%
  t() %>%
  scale() %>%
  t()
Heatmap(ht_sum,cluster_columns = F)
