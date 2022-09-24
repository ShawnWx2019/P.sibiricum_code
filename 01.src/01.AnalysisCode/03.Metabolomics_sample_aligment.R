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

load("../Data_cleaning/NEG/cleandata/object_neg.ms2.rds")
load("../Data_cleaning/POS/cleandata/object_pos.ms2.rds")

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
sample_info =
  object_origin %>%
  extract_sample_info() %>%
  dplyr::filter(sample_id != "Sample4_2" & sample_id != "Sample7_1")

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



plt = hclust(dist(t(expMat)))

plot(plt)

cor_trans = cor(expMat)

corrplot::corrplot(cor_trans,method = "number",hclust.method = 'complete',order = 'hclust')



# cor Transcript --------------------------------------------------------------------



expMat_trans <- read.table("../../../02.Data/02.Pacbio_raw/merged.fpkm.xls",header = T,sep = "\t")

expMat_trans <-
  expMat_trans %>%
  dplyr::select(!contains("9")) %>%
  dplyr::select(!all_of(c("HJ.5y.1","HJ.5y.2","HJ.7y.1"))) %>%
  dplyr::filter(rowSums(. > 1) > 1) %>%
  column_to_rownames("gene_id")

meta <- data.frame(
  row.names = colnames(expMat_trans),
  group = gsub("..$","",colnames(expMat_trans))
)

obj.pca2 <- pca(expMat_trans %>% as.matrix(),metadata = meta,center = T,scale = T,removeVar = .1)
x.name = data.frame(
  group = gsub("..$","",rownames(obj.pca2$rotated)),
  row.names  = rownames(obj.pca2$rotated)
)
plot_ly(x = obj.pca2$rotated$PC1,y = obj.pca2$rotated$PC2,z = obj.pca2$rotated$PC3,color = rownames(obj.pca2$rotated))

transcript <-
biplot(obj.pca2,colby = "group")

transcript+metabo

dist.trans = dist(t(expMat_trans))

plt = hclust(dist.trans)
plot(plt)


cor_trans = cor(expMat_trans)

corrplot::corrplot(cor_trans,method = "number",hclust.method = 'complete',order = 'hclust')

