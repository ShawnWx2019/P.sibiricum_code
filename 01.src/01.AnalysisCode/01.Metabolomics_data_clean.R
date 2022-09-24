######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: Data cleaning for metabolomics.
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################
getwd()

dir.create("../03.Progress/01.Metabolomics/Data_cleaning",showWarnings = F,recursive = T)

setwd("../03.Progress/01.Metabolomics/Data_cleaning")

library(tidymass)
library(tidyverse)
library(IMOtoolkits)
library(patchwork)
pacman::p_load(affy, parallel,ranger,caret,pcaMethods,ggplot2,tidyr,graphics,grDevices,Hmisc,gtools,cowplot,RColorBrewer,readr,plotly,stringr,GGally,dplyr,e1071,officer,bootstrap,pdist,metabolomics)

# functions for analysis --------------------------------------------------

## ggplot theme
theme1 =   theme(
  panel.border = element_rect(size = 1.5),
  axis.title = element_text(size = 14,color = 'black'),
  axis.text = element_text(size = 12,color = 'black')
)

# 01. NEG
# 1.1 batch effect detection --------------------------------------------------
load("../../../02.Data/01.Metabo_Raw/object.neg")
#> check positive

object

object %>% extract_sample_info()
sample_info_out <- read.csv("../../../02.Data/01.Metabo_Raw/Sample_information.csv")
#> correct the injection order
object <-
  object %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::select(sample_id) %>%
  dplyr::left_join(.,sample_info_out ,"sample_id") %>%
  arrange(injection.order)
object %>% extract_sample_info()
dir.create("NEG/raw",showWarnings = F,recursive = T)

save(object,file = "NEG/raw/object.rds")

#> batch effect
#>

plt_batch =
  object %>%
  extract_expression_data() %>%
  dplyr::select(contains("QC")) %>%
  pivot_longer(contains("QC"),values_to = "value",names_to = "QC_samples") %>%
  drop_na() %>%
  mutate(value = log2(value)) %>%
  ggplot(data = .,mapping = aes(x = QC_samples,y = value,color = QC_samples))+
  geom_boxplot()+
  ylab("log2(Raw peak area)")+
  xlab("")+
  theme_bw()+
  theme1+
  theme(
    axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1)
  )
# > have batch effect.
ggsave(plt_batch,filename = "NEG/Batch_detected.png",width = 8,height = 5,dpi = 300)


# mv and outlier detection ------------------------------------------------

plt_peak_dis =
  object %>%
  `+`(1) %>%
  log(10) %>%
  show_mz_rt_plot() +
  scale_size_continuous(range = c(0.01, 2))+
  theme1
ggsave(plt_peak_dis,filename = "NEG/peakdistribution.png",width = 13,height = 6,dpi = 300)


#> missing value

mv.plt <-
massdataset::show_sample_missing_values(object = object,
                                        color_by = "class",
                                        order_by = "injection.order",
                                        percentage = TRUE) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

ggsave(mv.plt,filename = "NEG/Missing_value.png",width = 15,height = 6,dpi = 300)

#> Sample_7_3 and Sample_9_3 ?
#>

outlier_samples =
  object %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

# add batch effect
object <-
  object %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::mutate(batch = rep(c("1","2"),c(7,23-6)))
object %>% extract_sample_info()
massqc_report(object = object,"NEG/data_quality_before_cleaning")

mv.plt2 <-
  massdataset::show_sample_missing_values(object = object,
                                          color_by = "class",
                                          order_by = "injection.order",
                                          percentage = TRUE) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

ggsave(mv.plt2,filename = "NEG/Missing_value_rm_outlier.png",width = 15,height = 6,dpi = 300)

massqc::massqc_report(object = object,
                      path = "NEG/data_quality_before_data_cleaning")

# remove noise ------------------------------------------------------------

qc_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)
object =
  object %>%
  mutate_variable_na_freq(according_to_samples = qc_id)
object <-
  object %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2)

#> outlier detective
massdataset::show_sample_missing_values(object = object,
                                        color_by = "class",
                                        order_by = "injection.order",
                                        percentage = TRUE) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

outlier_samples =
  object %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

outlier_table <-
  extract_outlier_table(outlier_samples)
## remove outliers
mv_out_name <-
outlier_table %>%
  apply(1, function(x){
    sum(x)
  }) %>%
  `>`(0) %>%
  which() %>% names()
object_mv_outlier <-
  object %>%
  activate_mass_dataset("expression_data") %>%
  dplyr::select(!all_of(mv_out_name))

## missing value imputation
object_impute <-
  impute_mv(object = object_mv_outlier, method = "knn")

norm_before <-
  IMO_plt_pca(obj = object_impute,tag = "group")


# normalization -----------------------------------------------------------

#> serrf
# set.seed(20164)
# object_SERRF <-
#   run_serrf(obj = object_impute,QC_tag = "QC",cluster_num = 8)
#
# norm_serrf <-
#   object_SERRF %>%
#   IMO_plt_pca(tag = "group")
#> svr
object_svr <-
  normalize_data(object_impute, method = "svr")

#> align batch
object_integration <-
  integrate_data(object_svr, method = "qc_mean")

intergrate_svr =
  object_integration %>%
  IMO_plt_pca(tag = 'group',interactive = F)


PCA.plt <- norm_before$plot + intergrate_svr$plot

ggsave("NEG/PCA_result.png",width = 16,height = 10,dpi = 300)

object.neg <- object_integration
dir.create("NEG/cleandata",showWarnings = F,recursive = T)
save(object_impute,file = "NEG/cleandata/object_impute.rds")
save(object_SERRF,file = "NEG/cleandata/object_SERRF.rds")
save(object.neg,file = "NEG/cleandata/object_neg.rds")
load("NEG/cleandata/object_neg.rds")
object.neg
# 2.0 POS -----------------------------------------------------------------
load("../../../02.Data/01.Metabo_Raw/object.pos")
#> check positive

object

object %>% extract_sample_info()
sample_info_out <- read.csv("../../../02.Data/01.Metabo_Raw/Sample_information.csv")
#> correct the injection order
object <-
  object %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::select(sample_id) %>%
  dplyr::left_join(.,sample_info_out,"sample_id") %>%
  arrange(injection.order) %>%
  dplyr::mutate(injection.order = 1:ncol(.))

dir.create("POS/raw",showWarnings = F,recursive = T)

save(object,file = "POS/raw/object.rds")

#> batch effect
#>

plt_batch =
  object %>%
  extract_expression_data() %>%
  dplyr::select(contains("QC")) %>%
  pivot_longer(contains("QC"),values_to = "value",names_to = "QC_samples") %>%
  drop_na() %>%
  mutate(value = log2(value)) %>%
  ggplot(data = .,mapping = aes(x = QC_samples,y = value,color = QC_samples))+
  geom_boxplot()+
  ylab("log2(Raw peak area)")+
  xlab("")+
  theme_bw()+
  theme1+
  theme(
    axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1)
  )

ggsave(plt_batch,filename = "POS/Batch_detected.png",width = 8,height = 5,dpi = 300)


# mv and outlier detection ------------------------------------------------

plt_peak_dis =
  object %>%
  `+`(1) %>%
  log(10) %>%
  show_mz_rt_plot() +
  scale_size_continuous(range = c(0.01, 2))+
  theme1
ggsave(plt_peak_dis,filename = "POS/peakdistribution.png",width = 13,height = 6,dpi = 300)


#> missing value

mv.plt <-
  massdataset::show_sample_missing_values(object = object,
                                          color_by = "class",
                                          order_by = "injection.order",
                                          percentage = TRUE) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

ggsave(mv.plt,filename = "POS/Missing_value.png",width = 15,height = 6,dpi = 300)

#> Sample_7_3 and Sample_9_3 ?
#>

outlier_samples =
  object %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

# remove outlier_samples
outlier_table <-
  extract_outlier_table(outlier_samples)

outlier_name <-
  outlier_table %>%
  apply(1, function(x){
    sum(x)
  }) %>%
  `>`(0) %>%
  which() %>%
  names()

# add batch effect

massqc_report(object = object,"POS/data_quality_before_cleaning")

mv.plt2 <-
  massdataset::show_sample_missing_values(object = object,
                                          color_by = "class",
                                          order_by = "injection.order",
                                          percentage = TRUE) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

ggsave(mv.plt2,filename = "POS/Missing_value_rm_outlier.png",width = 15,height = 6,dpi = 300)

massqc::massqc_report(object = object,
                      path = "POS/data_quality_before_data_cleaning")

# remove noise ------------------------------------------------------------

qc_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)
object =
  object %>%
  mutate_variable_na_freq(according_to_samples = qc_id)
object <-
  object %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2)

#> outlier detective
massdataset::show_sample_missing_values(object = object,
                                        color_by = "class",
                                        order_by = "injection.order",
                                        percentage = TRUE) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

outlier_samples =
  object %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

outlier_table <-
  extract_outlier_table(outlier_samples)
## remove outliers
mv_out_name <-
  outlier_table %>%
  apply(1, function(x){
    sum(x)
  }) %>%
  `>`(0) %>%
  which() %>% names()
object_mv_outlier <-
  object %>%
  activate_mass_dataset("expression_data") %>%
  dplyr::select(!all_of(mv_out_name))

## missing value imputation
object_impute <-
  impute_mv(object = object_mv_outlier, method = "knn")

norm_before <-
  IMO_plt_pca(obj = object_impute,tag = "group")


# normalization -----------------------------------------------------------

#> serrf
set.seed(20164)
object_SERRF <-
  run_serrf(obj = object_impute,QC_tag = "QC",cluster_num = 8)

norm_serrf <-
  object_SERRF %>%
  IMO_plt_pca(tag = "group")
#> svr
object_svr <-
  normalize_data(object_impute, method = "svr")

#> align batch
object_integration <-
  integrate_data(object_svr, method = "qc_mean")

intergrate_svr =
  object_integration %>%
  IMO_plt_pca(tag = 'group')


PCA.plt <- norm_before$plot + intergrate_svr$plot
object.pos <- object_integration
ggsave("POS/PCA_result.png",width = 16,height = 10,dpi = 300)
dir.create("POS/cleandata",showWarnings = F,recursive = T)
save(object_impute,file = "POS/cleandata/object_impute.rds")
save(object_SERRF,file = "POS/cleandata/object_SERRF.rds")
save(object.pos,file = "POS/cleandata/object_pos.rds")
load("POS/cleandata/object_pos.rds")
object.pos
object_origin <-
  merge_mass_dataset(
    x = object.neg,
    y = object.pos,
    sample_direction = "inner",
    variable_direction = "full",
    sample_by = "sample_id",
    variable_by = c("variable_id","mz","rt")
  ) %>%
  activate_mass_dataset("sample_info") %>%
  filter(class != "QC")
save(object_origin,file = "../Data_cleaning/object_origin_merge.rds")
