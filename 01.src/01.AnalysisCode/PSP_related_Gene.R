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


psp_gene_anno <- readxl::read_xlsx("../PSP_related_Gene.xlsx",sheet = 2)
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

c1_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`1`,]

c2_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`2`,]

psp_new_mat <- rbind(c6_mat,c1_mat,c2_mat)
# sum ---------------------------------------------------------------------


psp_mat_sum =
  psp_new_mat %>%
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
psp_mat_mean =
  psp_new_mat %>%
  select(-New_ID) %>%
  pivot_longer(
    !category,names_to = 'sample',values_to = 'value'
  ) %>%
  mutate(
    sample = gsub("..$","",sample)
  ) %>%
  group_by(category,sample) %>%
  summarise(mean = mean(value)) %>%
  pivot_wider(
    names_from = sample,values_from = mean
  )

ht_mean <-
  psp_mat_mean %>%
  column_to_rownames("category") %>%
  t() %>%
  scale() %>%
  t()
Heatmap(ht_mean,cluster_columns = F)


new <- read.delim("/Users/shawn/SynologyDrive/Project/11.HJ/01.Data/report/result/Annotation/HJ/new.xls",header = T,sep = "\t")


new_mat <- inner_join(
  data_frame(
    gene_id = new$gene_id,
    category = new$category,
  ),raw_fpkm
)

x_ht_mat = new_mat %>%
  column_to_rownames("gene_id") %>%
  filter(rowSums(.[-1] > 1) >2) %>%
  mutate(
    NewID = paste0("Gene",c(1:nrow(.)))
  )
ID_check = data.frame(
  old = rownames(x_ht_mat),
  new = x_ht_mat$NewID
)
ht_mat <-
  x_ht_mat%>%
  select(-category) %>%
  rownames_to_column("oldID") %>%
  column_to_rownames("NewID") %>%
  select(-oldID) %>%
  t() %>%
  scale() %>%
  t()

l_anno <- rowAnnotation(
  category = data.frame(
    row.names = x_ht_mat$NewID,
    cate = x_ht_mat$category
  ) %>% as.matrix(),
  col = list(
    category = c("galU" = "red","GMMP" = "green","PMM" = "purple","UXE" = "yellow")
  )
)

Heatmap(ht_mat,cluster_columns = F,row_km = 4,row_km_repeats = 2000,left_annotation = l_anno)

c6_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`6`,]

c1_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`1`,]

c2_mat <- psp_mat_clean[row_order(ht_mat_plt2)$`2`,]

psp_new_mat <- rbind(c6_mat,c1_mat,c2_mat)


head(ID_check)
