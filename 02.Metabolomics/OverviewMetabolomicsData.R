######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: overview
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################
dir.create("../03.Progress/01.Metabolomics/Overview")
setwd("../03.Progress/01.Metabolomics/Overview")

library(tidymass)
library(tidyverse)
library(IMOtoolkits)
library(PCAtools)
library(plotly)
# all feature PCA ---------------------------------------------------------

##> load massdataset
load("../Data_cleaning/NEG/cleandata/object_neg.ms2.rds")
load("../Data_cleaning/POS/cleandata/object_pos.ms2.rds")

##> merge neg and pos data
object_neg.ms2 <-
  object_neg.ms2 %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::mutate(sample_id = gsub("_neg","",sample_id))

object_pos.ms2 <-
  object_pos.ms2 %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::mutate(sample_id = gsub("_pos","",sample_id))

object_origin <-
  merge_mass_dataset(
    x = object_neg.ms2,
    y = object_pos.ms2,
    sample_direction = "inner",
    variable_direction = "full",
    sample_by = "sample_id",
    variable_by = c("variable_id","mz","rt")
  )
object_origin <-
  object_origin %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::filter(!class == "QC")
save(object_origin,file = "../Data_cleaning/object_origin_merge.rds")
##> extract data and pca analysis
sample_info =
  object_origin %>%
  extract_sample_info()
expMat <-
  object_origin %>%
  extract_expression_data()
write.csv(expMat,file = "expMat.csv")
meta <- data.frame(
  row.names = rownames(obj.pca$rotated),
  group = paste0("size-",c(3,3,3,5,1,4,2,2,6,1,"x",4,2,1,6,5))
)
expMat <-
  expMat %>% select(sample_info$sample_id) %>%
  mutate_if(is.numeric,.funs = function(x) log10(x+1)) %>% drop_na() %>% as.matrix()
obj.pca <- pca(expMat,metadata = meta,center = T,scale = T,removeVar = .1)

metabo <-
  biplot(obj.pca,colby = "group")

plot_ly(x = obj.pca$rotated$PC1,y = obj.pca$rotated$PC2,z = obj.pca$rotated$PC3,color = meta$group)

#> hclust analysis
plt = hclust(dist(t(expMat)))

plot(plt)

#> correlation analysis
cor_trans = cor(expMat)

corrplot::corrplot(cor_trans,method = "number",hclust.method = 'complete',order = 'hclust')
##> data vis.

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

##> data vis2 figS1
## load object
library(patchwork)
library(tidymass)
load("../Data_cleaning/NEG/raw/object.rds")
object_n_impute <- object
load("../Data_cleaning/POS/raw/object.rds")
object_p_impute <- object
load("../Data_cleaning/NEG/cleandata/object_neg.rds")
load("../Data_cleaning/POS/cleandata/object_pos.rds")
box_plt_batch = function(x,title) {
  x.order = x %>%
    extract_sample_info() %>%
    filter(group == "QC") %>%
    select(sample_id,injection.order) %>%
    arrange(injection.order) %>%
    pull(sample_id)
  box_plt <-
    x %>%
    extract_expression_data() %>%
    select(contains("QC")) %>%
    pivot_longer(contains("QC"),names_to = "sample",values_to = "peak") %>%
    mutate(peak = log10(peak+1),
           sample = factor(sample,levels = x.order)) %>%
    ggplot(.,aes(x = sample,y = peak,fill = sample)) + geom_boxplot(alpha = 0.6) +
    xlab("")+
    ylab("log10(peak area)")+
    ggtitle(title)+
    theme_bw()+
    theme(
      rect = element_rect(size = 1.5,color = 'black'),
      axis.title.y = element_text(size = 14,color = "black"),
      axis.text.x = element_text(size = 12,color = "black",angle = 90,vjust = .5,hjust = 1),
      axis.text.y = element_text(size = 12,color = "black"),
      title = element_text(size  = 15,colour = "black"),
      legend.position = "none"
    )
  return(box_plt)
}
n_before =
  object_n_impute %>% box_plt_batch(.,title = "negative before")
p_before =
  object_p_impute %>% box_plt_batch(.,title = "positive before")
n_after =
  object.neg %>% box_plt_batch(.,title = "negative after")
p_after =
  object.pos %>% box_plt_batch(.,title = "positive after")

figureS1 <-(p_before + p_after)/(n_before + n_after) + plot_annotation(tag_levels = "A")
ggsave(figureS1,filename = "../../../04.Result/02.Figures/FigS1.batch_effect_boxplot.png",width = 14,height = 7,dpi = 300)
ggsave(figureS1,filename = "../../../05.Report/FigS1.batch_effect_boxplot.png",width = 14,height = 7,dpi = 300)


library(PCAtools)

pos <- IMOtoolkits::IMO_plt_pca(obj = object.neg %>% activate_mass_dataset("expression_data") %>% drop_na(),tag = "group")

neg <- IMOtoolkits::IMO_plt_pca(obj = object.pos %>% activate_mass_dataset("expression_data") %>% drop_na(),tag = "group")

IMO_plt_rsd = function(obj_old,obj_new,QC_tag,x_loc = c(120,110),y_loc = c(15,110),interactive = F,showImage=T){
  # msg_yes = green$bold$italic
  # msg_no = red$bold$italic
  # msg_warning = yellow$bold$italic
  if (class(obj_old) != 'mass_dataset' | class(obj_new) != 'mass_dataset') {
    message(msg_no("error: Only obj which generated by massdataset accepted! Please check your input." ))
    return()
  }
  # calculate rsd for old object and new object.
  rsd_before <-
    obj_old %>%
    activate_mass_dataset("sample_info") %>%
    dplyr::filter(class == QC_tag) %>%
    extract_expression_data() %>%
    rownames_to_column("ID") %>%
    pivot_longer(contains(QC_tag),names_to = "tag",values_to = "value") %>%
    select(-tag) %>%
    group_by(ID) %>%
    summarise(
      raw.rsd = (sd(value,na.rm = T)/mean(value,na.rm = T))*100
    )
  rsd_after <-
    obj_new %>%
    activate_mass_dataset("sample_info") %>%
    dplyr::filter(class == QC_tag) %>%
    extract_expression_data() %>%
    rownames_to_column("ID") %>%
    pivot_longer(contains(QC_tag),names_to = "tag",values_to = "value") %>%
    select(-tag) %>%
    group_by(ID) %>%
    summarise(
      norm.rsd = (sd(value,na.rm = T)/mean(value,na.rm = T))*100
    )
  join_rsd <- inner_join(rsd_before,rsd_after,by = "ID") %>%
    mutate(
      norm_tag = case_when(
        norm.rsd <= 30 ~ "RSD ≤ 30",
        TRUE ~ "RSD > 30 "
      ),
      raw_tag = case_when(
        raw.rsd <= 30 ~ "RSD ≤ 30",
        TRUE ~ "RSD > 30 "
      )
    )
  #
  n1 = join_rsd %>% group_by(norm_tag) %>% summarise(num = n())

  rsd_plt = ggplot(data = join_rsd,mapping = aes(x = raw.rsd,y = norm.rsd,color = norm_tag)) +
    geom_point(size = 1.2,alpha =.8)+
    scale_color_manual(values = c("RSD ≤ 30" = "salmon","RSD > 30" = "grey30"))+
    geom_vline(xintercept = 30,linetype = "dashed",color = "red")+
    geom_hline(yintercept = 30,linetype = "dashed",color = "red")+
    geom_abline(slope = 1,linetype = "dashed",color = "red")+
    geom_label(
      data = data.frame(
        x = x_loc,
        y = y_loc,
        label = paste0("n=",c(n1[2,2],n1[1,2])),
        color = (c("RSD ≤ 30","RSD > 30"))
      ),
      mapping = aes(x = x,y = y,label = label,color= color),alpha = .8
    )+
    xlim(0,150)+
    ylim(0,150)+
    labs(x = "RSD of raw peak area",y = "RSD of svr normalized peak area")+
    theme_bw()+
    theme(
      line = element_line(size = 1.5,color = "black"),
      rect = element_rect(size = 1.5, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14,color = "black")
    )
  if (isTRUE(showImage)) {
    print(rsd_plt)
  }
  out = list(
    rsd_tbl = join_rsd,
    plot = rsd_plt
  )
  return(out)
}

rsd_neg <- IMO_plt_rsd(
  obj_old = object_n_impute,
  obj_new = object.neg,QC_tag = "QC"
)

rsd_pos <- IMO_plt_rsd(
  obj_old = object_p_impute,
  obj_new = object.pos,QC_tag = "QC"
)



figureS1 <-
  (pos$plot+ neg$plot)+ plot_annotation(tag_levels = "A")+ plot_layout(guides = "collect")&
  theme(legend.position =  'bottom',
        legend.title = element_blank())
figureS3 <-
  (rsd_pos$plot + rsd_neg$plot) + plot_annotation(tag_levels = "A")+ plot_layout(guides = "collect")&
  theme(legend.position =  'bottom',
        legend.title = element_blank())

ggsave(figureS3,filename = "../../../04.Result/02.Figures/FigS2.raw_data_pca.pdf",width = 14,height = 7.5,dpi = 300)
ggsave(figureS3,filename = "../05.Report/FigS2.raw_data_pca.pdf",width = 14,height = 8,dpi = 300)
