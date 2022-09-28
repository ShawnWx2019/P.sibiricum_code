######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: Metabolomics overview
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################

#dir.create("../03.Progress/01.Metabolomics/Overview/",showWarnings = F,recursive = T)

library(tidymass)
library(tidyverse)
library(PCAtools)
library(ComplexHeatmap)
library(scatterplot3d)
if(getwd() == "/Users/shawn/My_Repo/P.sibiricum/01.src") {
  setwd("../03.Progress/01.Metabolomics/Overview/")
} else {
  setwd("/Users/shawn/My_Repo/P.sibiricum/03.Progress/01.Metabolomics/Overview/")
}

# 01.Database ------------------------------------------------------------------
feature_anno <- readxl::read_xlsx("../Data_annotation_filter/Data_out/feature_anno_final.xlsx")

load("../../01.Metabolomics/Data_annotation_filter/Data_out/object_anno_merge.rds")
cate <- read.csv("../../01.Metabolomics/Data_annotation_filter/category.csv")

# Heatmap ---------------------------------------------------------------------


ExpMat <-
  left_join(
  data.frame(a = feature_anno$variable_id),
  object_anno_merge@expression_data %>% rownames_to_column("a")
) %>%
  distinct() %>%
  column_to_rownames("a")

ht_expmap_df <-
  ExpMat %>%
  rownames_to_column("metaID") %>%
  pivot_longer(!metaID,names_to = "sample",values_to = 'value') %>%
  mutate(sample = gsub("..$","",sample)) %>%
  group_by(metaID,sample) %>%
  summarise(mean = mean(value)) %>%
  pivot_wider(names_from = sample,values_from = mean) %>%
  column_to_rownames("metaID") %>%
  log10() %>%
  t() %>%
  scale() %>%
  t()

expmat_mean =
  ExpMat %>%
  rownames_to_column("metaID") %>%
  pivot_longer(!metaID,names_to = "sample",values_to = 'value') %>%
  mutate(sample = gsub("..$","",sample)) %>%
  group_by(metaID,sample) %>%
  summarise(mean = mean(value)) %>%
  pivot_wider(names_from = sample,values_from = mean) %>%
  column_to_rownames("metaID")

ht_expmat <-
  Heatmap(
  matrix = ht_expmap_df,
  cluster_columns = F,
  show_row_names = F,
  row_km = 8,
  row_km_repeats = 5000,
  name = " ",border = T,
  )
ht_expmat2 = draw(ht_expmat)
xxxx = row_order(ht_expmat2)
names(xxxx)

cluster_mat = map2_dfr(.x = xxxx,.y = names(xxxx),.f = function(.x,.y) {
  c_mat = ht_expmap_df[.x,] %>% as.data.frame %>% mutate(cluster = .y) %>% relocate(cluster,.before = Size_1)
  return(c_mat)
})

dot_tbl <-
cluster_mat %>%
  rownames_to_column('meta_ID') %>%
  pivot_longer(contains("Size"),names_to = "sample",values_to = "value") %>%
  mutate(cluster = factor(cluster,levels = names(xxxx))) %>%
  group_by(cluster,sample) %>%
  summarise(mean = mean(value)) %>%
  ungroup() %>%
  as.data.frame()

exp_table <- map2_dfr(.x = xxxx,.y = names(xxxx),.f = function(.x,.y) {
  c_mat = expmat_mean[.x,] %>% as.data.frame %>% mutate(cluster = .y) %>% relocate(cluster,.before = Size_1)
  return(c_mat)
})

exp_table_out <-
  left_join(exp_table %>% rownames_to_column("variable_id"),feature_anno) %>%
  left_join(.,cate %>% select(variable_id,Superclass,Class,Subclass))

writexl::write_xlsx(exp_table_out,"KM_annotation.xlsx")

View(c5)
trend <-
ggplot(dot_tbl,aes(x = sample,y = mean,color = cluster)) + geom_line(aes(group = cluster),size = 1.5) + facet_wrap(~cluster) +
  theme_bw()+
  theme(axis.text.x  = element_text(size = 12,color = "black",angle = 90),rect = element_rect(size = 1.5),axis.text.y =  element_text(color="black",size = 12) ,legend.position = "")

pdf("heatmap_kmeans.pdf",width = 6.5,height = 14)
ht_expmat2
dev.off()

pdf("kmeans_trend.pdf",width = 12,height = 9)
trend
dev.off()


# PCA ---------------------------------------------------------------------
load("../Data_cleaning/object_origin_merge.rds")
ExpMat1 <-
  object_origin %>%
  extract_expression_data()
ExpMat1 <-
  ExpMat1 %>%
  select(order(colnames(.))) %>% as.matrix()
metadata = data.frame(
  row.names = colnames(ExpMat1),
  Size = gsub("..$","",colnames(ExpMat1))
)
obj.pca = pca(
  mat = ExpMat1 %>% log10,metadata = metadata,center = T,scale = T,removeVar = .1
)

biplot(obj.pca,colby = "Size")

IMOtoolkits::IMO_plt_pca(obj = object_origin,tag = "group",center = T,scale = T,removeVar = .1,interactive = T)

pdf("pca3D.pdf",height = 7,width = 7)
scatterplot3d(
  obj.pca$rotated[,c(3,2,1)],
  xlab = paste0("PC3 (",round(obj.pca$variance[3],2),"%)"),
  ylab = paste0("PC2 (",round(obj.pca$variance[2],2),"%)"),
  zlab = paste0("PC1 (",round(obj.pca$variance[1],2),"%)"),
  pch = 16,angle = 30,
  box = T,cex.symbols = 2,lty.hide = 2,lty.grid = 2,
  type = 'p',color = rep(rainbow(6),c(3,3,3,2,2,2))
)
legend("topleft",paste0("Size",1:6),fill = rainbow(6))
dev.off()
