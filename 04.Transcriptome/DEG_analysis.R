######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: overview
#         Date: Aug 01, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################
library(tidyverse)
getwd()
dir.create("../03.Progress/02.Pacbio/Exp_data_split")

if(getwd() == "/Users/shawn/My_Repo/P.sibiricum/01.src") {
  setwd("../03.Progress/02.Pacbio/Exp_data_split")
} else {
  setwd("/Users/shawn/My_Repo/P.sibiricum/03.Progress/02.Pacbio/Exp_data_split")
}

library(DESeq2)
library(PCAtools)

# remove lncRNA -----------------------------------------------------------
FPKM <- read.table("../../../02.Data/02.Pacbio_raw/merged.fpkm.xls",header = T,sep = "\t")
count <- read.table("../../../02.Data/02.Pacbio_raw/merged.readcount.xls",header = T,sep = "\t")
transdecoder <- read.delim("../../../02.Data/02.Pacbio_raw/transdecoder.list.xls",sep = "\t",header = F) %>%
  mutate(
    V1 = str_extract(string = V1,pattern = "(?<=\\>).*\\d(?=\\:)")
  ) %>%
  setNames("gene_id")

expMat_clean <-
  inner_join(FPKM,transdecoder,"gene_id") %>%
  setNames(
    c("gene_id",
      "Size_1_1","Size_1_2","Size_1_3",
      "Size_2_1","Size_2_2","Size_2_3",
      "Size_3_1","Size_3_2","Size_3_3",
      "Size_4_1","Size_4_2",
      "Size_5_1","Size_5_2",
      "Size_6_1","Size_6_2")
  )
count_clean <-
  inner_join(count,transdecoder,"gene_id") %>%
  setNames(
    c("gene_id",
      "Size_1_1","Size_1_2","Size_1_3",
      "Size_2_1","Size_2_2","Size_2_3",
      "Size_3_1","Size_3_2","Size_3_3",
      "Size_4_1","Size_4_2",
      "Size_5_1","Size_5_2",
      "Size_6_1","Size_6_2")
  ) %>%
  mutate_if(is.numeric,ceiling)
write.table(count_clean,"../../../02.Data/02.Pacbio_raw/merged.readcount_clean.xls",row.names = F,sep = "\t")
write.table(expMat_clean,"../../../02.Data/02.Pacbio_raw/merged.FPKM_clean.xls",row.names = F,sep = "\t")


# PCA ---------------------------------------------------------------------
exp_mat <-
  expMat_clean %>%
  column_to_rownames('gene_id')

metadata <-
  data.frame(
    row.names = colnames(exp_mat),
    group = gsub("..$","",colnames(exp_mat)),
    stage = rep(c('Juvenal','Adult'),c(11,4))
  )
x_log <- exp_mat %>% +1 %>% log2()
mat.pca <- pca(mat = x_log,metadata = metadata,center = T,scale = T,removeVar = .1)

library(ggtree)
library(ade4)
dist = dist(exp_mat %>% t() %>% +1 %>% log2())
h_list = hclust(dist)

hc <- hclust(dist(mtcars))
hc
den = as.dendrogram(hc)
clus <- cutree(hc, 6)
g <- split(names(clus), clus)
library(ggtree)
library(ggplot2)
p <- ggtree(hc, linetype='dashed')
clades <- sapply(g, function(n) MRCA(p, n))

p <- groupClade(p, clades, group_name='subtree') + aes(color=subtree)

d <- data.frame(label = names(clus),
                cyl = mtcars[names(clus), "cyl"])

p %<+% d +
  layout_dendrogram() +
  geom_tippoint(aes(fill=factor(cyl), x=x+.5),
                size=5, shape=21, color='black') +
  geom_tiplab(aes(label=cyl), size=3, hjust=.5, color='black') +
  geom_tiplab(angle=90, hjust=1, offset=-10, show.legend=FALSE) +
  scale_color_brewer(palette='Set1', breaks=1:4) +
  theme_dendrogram(plot.margin=margin(6,6,80,6)) +
  theme(legend.position=c(.9, .6))



ggsave(tree.plt,filename = "ggtree.pdf",width = 9,height = 5.5)

write.tree(as.phylo(h_list),file = "hclust.nwk")

## visby ggtree

mytree <- as.dendrogram(h_list)

ggtree(mytree,ladderize = F) + geom_tiplab()

pdf("transcript_hclust.pdf",width = 8,height = 8)
plot(h_list)
dev.off()
#3D
library(scatterplot3d)
pdf("transcript_3D_PCA.pdf",width = 8,height = 8)
scatterplot3d(
  mat.pca$rotated[,c(1,2,3)],
  xlab = paste0("PC1 (",round(mat.pca$variance[1],2),"%)"),
  ylab = paste0("PC2 (",round(mat.pca$variance[2],2),"%)"),
  zlab = paste0("PC3 (",round(mat.pca$variance[3],2),"%)"),
  pch = 16,angle = 30,
  box = T,cex.symbols = 2,lty.hide = 2,lty.grid = 2,
  type = 'p',color = rep(rainbow(6),c(3,3,3,2,2,2))
)
legend("topleft",paste0("Size",1:6),fill = rainbow(6))
dev.off()

count_clean <-
  count_clean %>% column_to_rownames("gene_id") %>%
  as.matrix()

sample_info <-
  data.frame(
    row.names = colnames(count_clean),
    condition = rep(c('Juvenal','Adult'),c(11,4))
  )

dds <- DESeqDataSetFromMatrix(count_clean,sample_info,design = ~condition)
dds <- DESeq(dds)

res_all = results(dds,contrast = c("condition","Adult","Juvenal"))
DET =
  as.data.frame(res_all) %>%
  filter(
    abs(log2FoldChange) >= 0.26 & pvalue < 0.05 & !is.na(padj)
  )
writexl::write_xlsx(DET,"DET/AvsJ.xlsx")
# DEanalysis --------------------------------------------------------------

dir.create("DET",showWarnings = F,recursive = T)


sample_info <-
  data.frame(
    row.names = colnames(count_clean),
    condition = gsub("..$","",colnames(count_clean))
  )

dds <- DESeqDataSetFromMatrix(count_clean,sample_info,design = ~condition)
dds <- DESeq(dds)
left = rep(paste0("Size_",1:5),c(5:1 ))
right = rep(paste0("Size_",c(2:6,3:6,4:6,5,6,6)))
res.list = map2(.x = left,.y = right,.f = function(.x,.y) {
  res = results(dds,contrast = c("condition",.x,.y))
  DEG = subset(res,padj < 0.05 & abs(log2FoldChange) > 1) %>% as.data.frame() %>% rownames_to_column("Gene_ID")
  writexl::write_xlsx(x = DEG,path = paste0("DET/",.x,"_vs_",.y,".xlsx"))
  return(DEG)
})

## DEG_summary
DEG_summary = map2_dfr(.x = res.list,.y = paste0(left," vs. ",right),.f = function(.x,.y) {
  .x %>%
    mutate(
      regular = case_when(
        log2FoldChange > 0 ~ "up",
        log2FoldChange < 0 ~ "down"
      ),
      n = 1) %>%
    select(regular,n) %>%
    group_by(regular) %>%
    summarise(num = sum(n)) %>%
    mutate(group = gsub("Size_","Size",.y))

})

DEG_summary_plt <-
  ggplot(DEG_summary,aes(x = group,y = num,fill = regular)) +
  geom_bar(position="dodge",stat="identity",alpha = 0.65,color = "black")+
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()+
  theme(
    rect = element_rect(size = 1,color = "black"),
    axis.text.y = element_text(size = 12,color = "black",family = "DejaVu Sans"),
    axis.text.x = element_text(size = 12,color = "black",family = "DejaVu Sans",angle = 90,vjust = .5,hjust = 1),
    legend.title = element_blank(),
    axis.title = element_blank()
  )
ggsave(DEG_summary_plt,filename = "DET_summary.pdf",width = 9,height = 5.5,device = cairo_pdf)

## Get DEG_merge table

DET_merge_fpkm <- map_dfr(.x = res.list ,.f = function(.x) {
  .x %>% select(Gene_ID)
}) %>% unique() %>%
  inner_join(.,expMat_clean,c("Gene_ID" = "gene_id"))

DET_merge_count <- map_dfr(.x = res.list ,.f = function(.x) {
  .x %>% select(Gene_ID)
}) %>% unique() %>%
  inner_join(.,count_clean %>% as.data.frame() %>% rownames_to_column("gene_id"),c("Gene_ID" = "gene_id"))

cor_df <-
  DET_merge_fpkm %>%
  column_to_rownames("Gene_ID")
library(ComplexHeatmap)
cor_mat <- cor(x = cor_df %>% as.matrix(),method = "spearman")
col_fun = circlize::colorRamp2(breaks = c(0,1),colors = c("white","purple"))
pdf("corr_DETs.pdf",width = 7.2,height = 7)
Heatmap(cor_mat,
        cluster_rows = F,
        cluster_columns = F,
        col = col_fun,
        name = "r",border = T)
dev.off()



# mFuzz -------------------------------------------------------------------

mfuzz_input <-
  DET_merge_fpkm %>%
  pivot_longer(!Gene_ID,names_to = 'sample',values_to = 'value') %>%
  mutate(sample = gsub("..$","",sample)) %>%
  group_by(Gene_ID,sample) %>%
  summarise(mean = mean(value)) %>%
  pivot_wider(names_from = sample,values_from = mean) %>%
  column_to_rownames("Gene_ID")

obj_mfuzz <- new("ExpressionSet",exprs = mfuzz_input %>% as.matrix())
dt.s <- standardise(obj_mfuzz);
m1 <- mestimate(dt.s);
set.seed(20220921)
cl <- mfuzz(dt.s,c=12,m = m1)
png("mFuzz_Transcript.png",width = 4600,height = 1200,res = 300)
print({
  mfuzz.plot2(dt.s,cl,mfrow = c(2,6),x11 = FALSE,time.labels = 1:6,xlab = "",ylab = "Scaled peak")
})
dev.off()
mat_cluster <-
  cl$cluster %>%
  as.data.frame() %>%
  rownames_to_column("Gene_ID") %>%
  setNames(c("Gene_ID","cluster"))  %>%
  inner_join(mfuzz_input %>%
               as.data.frame %>%
               rownames_to_column("Gene_ID"),
             "Gene_ID") %>%
  arrange(cluster)
c2_id = mat_cluster %>% filter(cluster == 2) %>% pull(Gene_ID)

library(clusterProfiler)
GO_database <- readxl::read_xlsx("../../../02.Data/02.Pacbio_raw/GOdb.xlsx")
head(GO_database)
GO_db_hj <-
  GO_database %>%
  pivot_longer(!GeneID,names_to = "label",values_to = "value") %>%
  drop_na() %>%
  select(GeneID,value) %>%
  mutate(
    TERM = str_split(value,pattern = ": ",n = 2,simplify = T)[,1],
    NAME = str_split(value,pattern = ": ",n = 2,simplify = T)[,2]
  ) %>%
  select(GeneID,TERM,NAME)
t2g = GO_db_hj %>% select(TERM,GeneID)
t2n = GO_db_hj %>% select(TERM,NAME)

ego <- enricher(
  gene = c2_id,TERM2GENE = t2g,TERM2NAME = t2n,pvalueCutoff = 1,qvalueCutoff = 1
)
dotplot(ego)
write.table(t2g,file = "../../../02.Data/HJ.GO.txt",row.names = F,sep = "\t")
writexl::write_xlsx(mat_cluster,"mat_cluster.xlsx")
