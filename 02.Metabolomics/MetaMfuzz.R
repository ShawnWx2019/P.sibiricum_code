######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: mfuzz metabolomics
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################

dir.create("../03.Progress/01.Metabolomics/mFuzz")

if(getwd() == "/Users/shawn/My_Repo/P.sibiricum/01.src") {
  setwd("../03.Progress/01.Metabolomics/mFuzz")
} else {
  setwd("/Users/shawn/My_Repo/P.sibiricum/03.Progress/01.Metabolomics/mFuzz")
}

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require('Mfuzz')) BiocManager::install('Mfuzz');

library(Mfuzz)
library(edgeR)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
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



# DAM analysis ------------------------------------------------------------
library(MDAtoolkits)
library(tidyverse)
library(rstatix)
expmat2 <-
  expmat %>%
  distinct %>%
  column_to_rownames("name")
##> DAM analysis
tmpx <-
  MDAtoolkits::DM_analysis(
    x = expmat %>% distinct() %>% dplyr::rename("CompoundID" = "name") %>% dplyr::select(order(colnames(.))),
    right_index = c(1,2,3,4,5,6,7,8,9,10,11),
    left_index = c(12,13,14,15),
    right = "Juvenal",
    left = "Adult",
    method = "anova-test",
    method2 = "opls-da"
  )
DAM_tbl <-tmpx$DAM_tbl
a_vs_j_result <-
  DAM_tbl %>%
  filter(pvalue < 0.05  & abs(log2fc)  > 1) %>%
  left_join(annotation,c("CompoundID" = "variable_id")) %>%
  arrange(desc(abs(log2fc)))
a_vs_j_result$CompoundID %>% unique() %>% length()
writexl::write_xlsx(x = a_vs_j_result,"../../../04.Result/01.Tables/Adult_vs_Juvenal.xlsx")

expmat3 <- expmat %>%
  pivot_longer(!name,names_to = 'sample',values_to = 'peak') %>%
  mutate(sample = gsub("..$","",sample)) %>%
  group_by(name,sample) %>%
  summarise(mean = mean(peak)) %>%
  pivot_wider(names_from = sample,values_from = mean)

b_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("salmon","purple")),
                                             labels = c("Juvenal","Adult"),
                                             labels_gp = gpar(col = "white", fontsize = 10)))

a_vs_j_result2 <- readxl::read_xlsx("../../../04.Result/01.Tables/Adult_vs_Juvenal.xlsx")
##> heatmap
library(ComplexHeatmap)
library(circlize)
key_metabo <-
  x.order_lab %>%
  select(CompoundID,Compound.name) %>% distinct() %>%
  group_by(CompoundID) %>%
  top_n(1) %>%
  ungroup() %>%
  mutate(n = 1:nrow(.)) %>%
  inner_join(.,tmp_id,by = c("CompoundID" = "name")) %>%
  select(name2,n)

r_anno = rowAnnotation(
  foo = anno_mark(at = key_metabo %>% pull(n), labels = key_metabo %>% pull(name2)))


ht_cyc_path_2 =
  a_vs_j_result2 %>% select(CompoundID) %>%
  left_join(.,expmat3 %>% select(order(colnames(.))),c("CompoundID" = "name")) %>%
  distinct() %>%
  arrange(CompoundID) %>%
  column_to_rownames("CompoundID") %>%
  as.matrix() %>% log2() %>%
  Heatmap(
    .,
    name = "log2(peak area)",
    show_row_names = F,cluster_columns = F,
    col = colorRamp2(c(min(log2(expmat[,-1])), max(log2(expmat[,-1]))), c("yellow", "#FF1493")),
    border_gp = gpar(color = "black",size = 1),
    column_km  = 2,
    bottom_annotation = b_anno
  )
ht_cyc_path_3 = a_vs_j_result2 %>%
  select(CompoundID,log2fc) %>%
  distinct() %>%
  arrange(CompoundID) %>%
  column_to_rownames("CompoundID") %>%
  Heatmap(
    .,
    show_row_names = F,cluster_rows = F,
    name = "log2 (fold channge)",
    col = colorRamp2(c(-4,-2, 0, 2,4) ,c("blue","green", "white", "purple","red")),
    border_gp = gpar(color = "black",size = 1)
  )
ht_cyc_tag =
  ht_cyc_path_4 = a_vs_j_result2 %>%
  select(CompoundID,pvalue) %>%
  distinct() %>%
  arrange(CompoundID) %>%
  column_to_rownames("CompoundID") %>%
  mutate(
    pvalue = case_when(
      pvalue <= 0.0001 ~ "****",
      pvalue > 0.0001 & pvalue <= 0.001 ~ "***",
      pvalue > 0.001 & pvalue <= 0.01 ~ "**",
      pvalue > 0.01 & pvalue <= 0.05 ~ "*",
      TRUE ~ ""
    )
  )
ht_cyc_path_4 = a_vs_j_result2 %>%
  select(CompoundID,pvalue) %>%
  distinct() %>%
  arrange(CompoundID) %>%
  column_to_rownames("CompoundID") %>%
  Heatmap(
    .,
    show_row_names = F,cluster_rows = F,
    name = "pvalue",
    #row_order = eko_ht_order,
    col = colorRamp2(c(0.05,0), c("white", "salmon")),
    border_gp = gpar(color = "black",size = 1),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(ht_cyc_tag[i, j], x, y, gp = gpar(fontsize = 10))
    }
  )
ht_cyc_path_5 = a_vs_j_result2 %>%
  select(CompoundID,FDR) %>%
  distinct() %>%
  arrange(CompoundID) %>%
  column_to_rownames("CompoundID") %>%
  Heatmap(
    .,
    show_row_names = F,cluster_rows = F,
    name = "FDR",
    #row_order = eko_ht_order,
    col = colorRamp2(c(0.05,0), c("white", "#CC6699")),
    border_gp = gpar(color = "black",size = 1)
  )


ht_cyc_path_6 = a_vs_j_result2 %>%
  select(CompoundID,VIP) %>%
  distinct() %>%
  arrange(CompoundID) %>%
  column_to_rownames("CompoundID") %>%
  Heatmap(
    .,
    show_row_names = T,cluster_rows = F,
    name = "VIP",
    col = colorRamp2(c(2,0), c("cyan", "white")),
    border_gp = gpar(color = "black",size = 1)
  )
ht_cyc_enrich = ht_cyc_path_1+ht_cyc_path_2+ht_cyc_path_3+ht_cyc_path_4+ht_cyc_path_5+ht_cyc_path_6
pdf(file = paste0("../../../03.Progress/01.Metabolomics/mFuzz/DAM2.pdf"),width = 10.224,height = 16)
draw(ht_cyc_enrich)
invisible(dev.off())
png(file = paste0("../03.Progress/01.Metabolomics/mFuzz/DAM.png"),width = 2500,height = 2200,res = 300)
draw(ht_cyc_enrich)
invisible(dev.off())


# vocanoplot --------------------------------------------------------------

Vocanoplot = function(x,title){
  x %>%
    mutate(group = case_when(
      pvalue < pval_cut & VIP > VIP_cut & log2fc > log2fc_cut & FDR < qval_cut ~ "up",
      pvalue < pval_cut & VIP > VIP_cut & log2fc < -log2fc_cut & FDR < qval_cut~ "down",
      TRUE ~ "not significant"
    ),
    log10pval = -log10(pvalue),
    log10qval = -log10(FDR)
    )-> vocdata
  if(qval_cut <= pval_cut) {yinter = qval_cut} else {yinter = pval_cut}
  p = ggplot(data = vocdata,aes(x = log2fc,y = log10pval))+
    geom_point(aes(color = group,size = VIP),alpha = .7)+
    scale_color_manual(values = c("blue","grey","red"))+
    scale_size(range = c(0.05,4))+
    theme_bw()+
    ylab(expression(paste(-log10,("p-value"))))+
    xlab(expression(paste(log2,"(fold change)")))+
    xlim(-max(vocdata$log2fc)-1,max(vocdata$log2fc)+1)+
    ylim(0,max(vocdata$log10pval)+1)+
    ggtitle(title)+
    annotate("text",x = max(vocdata$log2fc)-1,y = max(vocdata$log10pval)-1,
             label = paste0("up regular:",nrow(filter(vocdata,group == "up"))),
             color = "red",size = 5)+
    annotate("text",x = -max(vocdata$log2fc)+1,y = max(vocdata$log10pval)-1,
             label = paste0("down regular:",nrow(filter(vocdata,group == "down"))),
             color = "blue",size = 5)+
    geom_hline(yintercept = -log10(yinter),lty=3,col="black",lwd=0.5)+
    geom_vline(xintercept = log2fc_cut,lty=3,col="black",lwd=0.5)+
    geom_vline(xintercept = -log2fc_cut,lty=3,col="black",lwd=0.5)
  theme(
    axis.title = element_text(size = 14,colour = 'black'),
    rect = element_rect(size = 1.5),
    axis.text = element_text(size = 12,colour = 'black'),
    title = element_text(size = 14,color = 'black'),
    text = element_text(size  = 12)
  )
  return(p)
}

pval_cut <- 0.05
VIP_cut <- 0
log2fc_cut <- 1
qval_cut <- 1
voc.plot <-
  Vocanoplot(DAM_tbl, "Adult vs Juvenal")
ggsave(voc.plot,filename = "../../../03.Progress/01.Metabolomics/mFuzz/vocanoplot.png",width = 7,height = 7)
ggsave(voc.plot,filename = "../../../03.Progress/01.Metabolomics/mFuzz/vocanoplot.pdf",width = 7,height = 7)

DAM_tbl


# barplot -----------------------------------------------------------------

a = readxl::read_xlsx("../../../02.Data/01.Metabo_Raw/key_metabolites.xlsx")
m_name = a$name2
key_mat <-
  inner_join(a,expmat3) %>%
  select(-name) %>%
  select(order(colnames(.))) %>%
  pivot_longer(!name2,names_to = "sample",values_to = 'value') %>%
  mutate(group = case_when(
    sample == "Size_5" | sample == "Size_6" ~ "Adult",
    TRUE ~ "Juvenal"
  ))

ggplot(
  data = key_mat,
  mapping = aes( x = sample,y = value,fill = group)
) + geom_bar(stat = 'identity')+facet_grid(~ name2)+
  theme_bw()+
  theme(
    axis.title = element_text(size = 14,color = "black"),
    axis.text.x = element_text(size = 12,colour = "black"),
    panel.border = element_rect(size = 1,colour = 'black')

  )


pltxx = map(.x = m_name,.f = function(.x) {
  ggplot(
    data = key_mat %>% filter(name2 == .x),
    mapping = aes( x = sample,y = value,fill = group)
  ) + geom_bar(stat = 'identity')+
    ggtitle(.x)+
    xlab("")+
    ylab("peak area")+
    theme_bw()+
    theme(
      axis.title = element_text(size = 14,color = "black"),
      axis.text.x = element_text(size = 12,colour = "black",angle = 90),
      axis.text.y = element_text(size = 12,colour = "black"),
      panel.border = element_rect(size = 1.5,colour = 'black'),
      title = element_text(size = 15,face = 'bold',color = 'black')
    )
})


library(patchwork)
xx_out <-
  pltxx[[1]]+pltxx[[2]]+pltxx[[3]]+pltxx[[4]]+pltxx[[5]]+pltxx[[6]]

ggsave(xx_out,filename = "key_metabolite.pdf",width = 12,height = 8)



# add classififcation -----------------------------------------------------

xx.exp <- readxl::read_xlsx("../../../04.Result/01.Tables/TableS4.DAM.xlsx",sheet = 2)
xx.classification <- readxl::read_xlsx("../../../04.Result/01.Tables/TableS3.Compound_classification.xlsx",sheet = 2)

xx.exp_out <-
  xx.exp %>%
  left_join(.,xx.classification,by = c( "CompoundID" = "variable_id"))
writexl::write_xlsx(xx.exp_out,path = "../../../04.Result/01.Tables/Adult_vs_Juvenal2.xlsx")



