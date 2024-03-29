###########################################################
#      Prj: Multi-omics data analysis of P.sibiricum.
#      Assignment: Phloretin
#      Date: Mar 21, 2022
#      Author: Shawn Wang <shawnwang2016@126.com>
#      Location: HENU, Kaifeng, Henan, China
###########################################################

if (getwd()!="/Users/shawn/My_Repo/P.sibiricum/01.src") {
  setwd("~/My_Repo/P.sibiricum/03.Progress/03.Transcriptome/Phloretin/")
} else {
  setwd("../03.Progress/03.Transcriptome/Phloretin")
}

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(rstatix)

# import data -------------------------------------------------------------

Phloretin_gene_anno <- readxl::read_xlsx("Phloretin.xlsx")
raw_fpkm <- readxl::read_xlsx("../FPKM.xlsx")

Phloretin_fpkm_raw <-
  inner_join(
    x = data.frame(
      gene_id = Phloretin_gene_anno$gene_id,
      category = Phloretin_gene_anno$category
    ),
    y = raw_fpkm,
    by = "gene_id"
  ) %>%
  filter(
    rowSums(.[,-c(1,2)] > 1) >2
  )
DET <- readxl::read_xlsx("../../02.Pacbio/Exp_data_split/DET/FPKM_merge_all.xlsx")
Phloretin_fpkm_raw1 = inner_join(Phloretin_fpkm_raw,
                                 data.frame(gene_id = DET$Gene_ID))
sig_result <-
  Phloretin_fpkm_raw1 %>%
  pivot_longer(contains('Size'),names_to = 'Size',values_to = 'value') %>%
  mutate(Size = gsub(pattern = "..$",replacement = "",Size)) %>%
  mutate(group = case_when(
    str_extract(Size,pattern = ".$") %>% as.numeric() <= 4 ~ "Juvenile",
    TRUE ~ "Adult"
  )) %>%
  group_by(category) %>%
  anova_test(value ~ group,detailed = T)
write.table(anno_result,"../../../05.Report/TableS12.phloretin.anno.xlsx")
sig_result2 <-
  Phloretin_fpkm_raw1 %>%
  pivot_longer(contains('Size'),names_to = 'Size',values_to = 'value') %>%
  mutate(Size = gsub(pattern = "..$",replacement = "",Size)) %>%
  mutate(group = case_when(
    str_extract(Size,pattern = ".$") %>% as.numeric() <= 4 ~ "Juvenile",
    TRUE ~ "Adult"
  )) %>%
  group_by(category) %>%
  wilcox_test(value ~ group,detailed = T)


# merge by sum values. ----------------------------------------------------

##> Merge all isoforms which annotated to the same 'gene'

Phloretin_exp_mat_long <-
  Phloretin_fpkm_raw1 %>%
  pivot_longer(contains('Size'),names_to = 'Size',values_to = 'value') %>%
  mutate(Size = gsub(pattern = "..$",replacement = "",Size)) %>%
  group_by(category,Size) %>%
  summarise(sum = sum(value))

##> For each Sizes of Juvenal ar Adult stage.
Phloretin_exp_mat_detail <-
  Phloretin_exp_mat_long %>%
  pivot_wider(names_from = Size,values_from = sum)
library(rstatix)
##> Gathered all Size and calculate the log2fc of A vs J
Phloretin_exp_mat_JvsA <-
  Phloretin_exp_mat_long %>%
  mutate(group = case_when(
    str_extract(Size,pattern = ".$") %>% as.numeric() <= 4 ~ "Juvenal",
    TRUE ~ "Adult"
  )) %>%
  group_by(category,group) %>%
  summarise(mean = mean(sum)) %>%
  pivot_wider(names_from = group,values_from = mean) %>%
  column_to_rownames("category") %>%
  mutate(log2fc = log2(Adult/Juvenal),
         FC = Adult/Juvenal)



col_fun = circlize::colorRamp2(c(-1,0,1),colors = c("blue","white","red"))

pdf("log2fc_Phloretin.pdf",width = 4,height = 9)

Heatmap(
  Phloretin_exp_mat_JvsA %>% select(log2fc),
  col = col_fun,
  name = "Adult vs Juvenal"
)
dev.off()

gene_anno <-
  inner_join(
    data.frame(
      gene_id = Phloretin_fpkm_raw$gene_id
    ),
    Phloretin_gene_anno
  ) %>% as_tibble() %>%
  left_join(
    .,data.frame(
      gene_id = Phloretin_fpkm_raw1$gene_id,
      `DET set` = TRUE
    )
  ) %>%
  mutate(
    DET.set =  case_when(
      is.na(DET.set) ~ FALSE,
      TRUE ~ TRUE
    )
  )
writexl::write_xlsx(gene_anno,"../../../05.Report/TableS11.Phloretin.xlsx")

