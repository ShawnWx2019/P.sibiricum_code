######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: mfuzz metabolomics
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################

dir.create("../03.Progress/01.Metabolomics/mFuzz")
setwd("../03.Progress/01.Metabolomics/mFuzz")
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require('Mfuzz')) BiocManager::install('Mfuzz');

library(Mfuzz)
library(edgeR)
library(tidyverse)
library(ComplexHeatmap)
category <- read.csv("../Data_annotation_filter/category.csv")
annotation <- readxl::read_xlsx("../Data_annotation_filter/Data_out/feature_anno_final.xlsx")
expmat <- read.csv("../Data_annotation_filter/Data_out/Peak_Table.csv") %>%
  select(-mz,-rt) %>%
  inner_join(annotation %>% select(variable_id),c("name" = "variable_id"))

expmat_mean <-
  expmat %>%
  pivot_longer(!name,names_to = "sample",values_to = "peak") %>%
  mutate(sample = gsub("..$","",sample)) %>%
  group_by(name,sample) %>%
  summarise(value = mean(peak)) %>%
  pivot_wider(names_from = sample,values_from = value) %>%
  column_to_rownames("name")

obj_mfuzz <- new("ExpressionSet",exprs = expmat_mean %>% as.matrix())
dt.s <- standardise(obj_mfuzz);
m1 <- mestimate(dt.s);
set.seed(20220901)
cl <- mfuzz(dt.s,c=15,m = m1)
print({
  mfuzz.plot2(dt.s,cl,mfrow = c(5,3),x11 = FALSE,time.labels = 1:6,xlab = "",ylab = "Scaled peak")
})
mat_cluster <-
  cl$cluster %>%
  as.data.frame() %>%
  rownames_to_column("Gene_ID") %>%
  setNames(c("Gene_ID","cluster"))  %>%
  inner_join(expmat_mean %>%
               as.data.frame %>%
               rownames_to_column("Gene_ID"),
             "Gene_ID") %>%
  arrange(cluster)
c7 <- mat_cluster %>% filter(cluster == 6)
plt_mt_7 <- c7 %>% column_to_rownames("Gene_ID") %>% select(-cluster) %>% log10() %>% t() %>% scale() %>% t()
Heatmap(plt_mt_7,show_row_names = F,cluster_columns = F)
c7 %>% left_join(category,c("Gene_ID" = "variable_id"))
