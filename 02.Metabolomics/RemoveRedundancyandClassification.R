######################################################################
#         Prj: Multi-omics data analysis of P.sibiricum.
#         Assignment: Remove redundancy and classification
#         Date: Mar 21, 2022
#         Author: Shawn Wang <shawnwang2016@126.com>
#         Location: HENU, Kaifeng, Henan, China
######################################################################
getwd()
library(tidymass)
library(tidyverse)
dir.create("../03.Progress/01.Metabolomics/Data_annotation_filter",showWarnings = F,recursive = T)
setwd("../03.Progress/01.Metabolomics/Data_annotation_filter")
load("../Data_cleaning/NEG/cleandata/object_neg_anno.rds")
load("../Data_cleaning/POS/cleandata/object_pos_anno.rds")

# change sample names for next steps.



# remove no annotation features -------------------------------------------
dir.create("high_confidence")
object_neg_clean_high <-
  object_neg_anno %>%
  activate_mass_dataset("annotation_table") %>%
  filter(!is.na(Level)) %>%
  filter(Level == 1 | Level == 2) %>%
  activate_mass_dataset("sample_info")

object_pos_clean_high <-
  object_pos_anno %>%
  activate_mass_dataset("annotation_table") %>%
  filter(!is.na(Level)) %>%
  filter(Level == 1 | Level ==2)%>%
  activate_mass_dataset("sample_info")

object_anno_merge <-
  merge_mass_dataset(
    x = object_neg_anno,
    y = object_pos_anno,
    sample_direction = "inner",
    variable_direction = "full",
    sample_by = "sample_id",
    variable_by = c("variable_id","mz","rt")
  ) %>%
  activate_mass_dataset("sample_info") %>%
  filter(class != "QC")
save(object_anno_merge,file = "Data_out/object_anno_merge.rds")
object_high <-
  merge_mass_dataset(
    x = object_neg_clean_high,
    y = object_pos_clean_high,
    sample_direction = "inner",
    variable_direction = "full",
    sample_by = "sample_id",
    variable_by = c("variable_id","mz","rt")
  ) %>%
  activate_mass_dataset("sample_info") %>%
  filter(class != "QC")
object_high <-
  object_high %>%
  activate_mass_dataset(what = "annotation_table") %>%
  group_by(Compound.name) %>%
  filter(Level == min(Level)) %>%
  filter(Total.score == max(Total.score)) %>%
  slice_head(n = 1)

save(object_high,file = "high_confidence/object_high.rda")
report_parameters(object_high,path = "high_confidence")
# with level 3
dir.create("all_features")
object_neg_clean <-
  object_neg_anno %>%
  activate_mass_dataset("annotation_table") %>%
  filter(!is.na(Level)) %>%
  activate_mass_dataset("sample_info")

object_pos_clean <-
  object_pos_anno %>%
  activate_mass_dataset("annotation_table") %>%
  filter(!is.na(Level)) %>%
  activate_mass_dataset("sample_info")

object_full <-
  merge_mass_dataset(
    x = object_neg_clean,
    y = object_pos_clean,
    sample_direction = "inner",
    variable_direction = "full",
    sample_by = "sample_id",
    variable_by = c("variable_id","mz","rt")
  ) %>%
  activate_mass_dataset("sample_info") %>%
  filter(class != "QC")
object_full<-
  object_full %>%
  activate_mass_dataset(what = "annotation_table") %>%
  group_by(Compound.name) %>%
  filter(Level == min(Level)) %>%
  filter(Total.score == max(Total.score)) %>%
  group_by(variable_id) %>%
  filter(Level == min(Level)) %>%
  filter(Total.score == max(Total.score)) %>%
  slice_head(n = 1)
object_full2 <-
  object_full %>%
  activate_mass_dataset("annotation_table") %>%
  mutate(
    judge = case_when(
      Level != 3 ~ 1,
      Level == 3 & Adduct == "(M-H)-" ~ 1,
      Level == 3 & Adduct == "(M+H)+" ~ 1,
      TRUE ~ 2
    )
  ) %>%
  filter(judge == 1)

save(object_full2,file = "all_features/object_full_feature.rda")
save(object_full,file = "all_features/object_full_feature2.rda")

report_parameters(object_full2,path = "all_features")


# output result -----------------------------------------------------------
dir.create("Data_out",showWarnings = F,recursive = T)
massdataset::export_mass_dataset(object = object_full,file_type = "csv",ms2_file_type = "mgf",path = "Data_out/")
massdataset::export_mass_dataset4metdna(object_full,path = "Data_out/")
dir.create("Data_out/NEG",showWarnings = F,recursive = T)
massdataset::export_mass_dataset4metdna(object = object_neg_anno,path = "Data_out/NEG/")
massdataset::export_mass_dataset4metdna(object = object_pos_anno,path = "Data_out/POS/")

full_anno_remove_reduandant <-
  object_full2 %>%
  extract_annotation_table()
writexl::write_xlsx(full_anno_remove_reduandant,path = "Data_out/feature_anno_remove_redundancy.xlsx")


# add inchikey ------------------------------------------------------------

full_anno_remove_reduandant <- readxl::read_xlsx("Data_out/feature_anno_remove_redundancy.xlsx")
load("~/SynologyDrive/database/02.MS/MSdb/labID2INCHIKEY.rda")
load("~/SynologyDrive/database/02.MS/MSdb/fiehn_hilic_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/hmdb_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/kegg_ms1_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/massbank_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/mona_database0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/snyder_database_hilic0.0.3.rda")
load("~/SynologyDrive/database/02.MS/MSdb/KNApSAcK_ms1_database.rda")
load("~/SynologyDrive/database/02.MS/MSdb/plantcyc_ms1_database0.0.2.rda")
load("~/SynologyDrive/database/02.MS/MSdb/RIKEN_PlaSMA_database0.0.1.rda")
load("~/SynologyDrive/database/02.MS/MSdb/Natural_products_database_v1.rds")

full_anno_inchi <-
  full_anno_remove_reduandant %>%
  left_join(lab2inchi) %>%
  left_join(object_full2@variable_info %>% select(variable_id,mz,rt),"variable_id")
view(full_anno_inchi)

##> classification
library(MDAtoolkits)
step2_result =
  data.frame(
    query = full_anno_inchi$Compound.name,
    cid = full_anno_inchi$PubChem.ID
  ) %>%
  filter(!is.na(cid))
step3_result <- mda_pubchem_crawler(
  cid_info = step2_result,type = "multiple",core_num = 8
)

final_anno <-
  full_anno_inchi %>%
  left_join(.,step3_result %>% select(CID,InChIKey) %>% mutate(CID = as.character(CID)),by = c("PubChem.ID" = "CID")) %>%
  mutate(INCHIKEY = case_when(
    is.na(INCHIKEY) ~ InChIKey,
    TRUE ~ INCHIKEY
  )) %>%
  select(-InChIKey)
writexl::write_xlsx(final_anno,path = "Data_out/feature_anno_final.xlsx")
tmp_x = readxl::read_xlsx("../../../04.Result/01.Tables/TableS1.Metabolomics_data.xlsx")
tmp_x =
  tmp_x %>%
  mutate(ion = case_when(
    str_detect(string = variable_id,"_POS") ~ "POS",
    str_detect(string = variable_id,"_NEG") ~ "NEG"
  ))
tmp_x %>% filter(ion == "POS") %>% select(variable_id,ion) %>% distinct() %>%  nrow()
tmp_x %>% filter(ion == "NEG") %>% select(variable_id,ion) %>% distinct() %>%  nrow()
source("~/My_Repo/MyBioScript/01.R/ShawnRToolkit.R")

a = BatchReadTable(path = "Data_out/classifire/",type = "csv",pattern = "classyfire.*",sep = ",",header = T )

classyfication <- bind_rows(a$file)

percent.tbl <- left_join(
  data.frame(
    variable_id = final_anno$variable_id,
    name = final_anno$Compound.name,
    InChIKey = final_anno$INCHIKEY
  ),classyfication
)
write.csv(percent.tbl,file = "category.csv",row.names = F)
mda_sum_categroy2 = function(cate_df,cate_type,mask_cutoff = 5,mask_cutoff2 = 3){
  categroy = cate_df %>%
    dplyr::rename("mda_class" = cate_type) %>%
    select(compound_id,mda_class) %>%
    mutate(
      n = 1
    ) %>%
    filter(!is.na(mda_class)) %>%
    group_by(mda_class) %>%
    summarise(count = sum(n)) %>%
    filter(mda_class != "")

  category_major = categroy %>%
    filter(count >= mask_cutoff)

  category_minor = categroy %>%
    filter(count < mask_cutoff) %>%
    mutate(
      class_less = map2_chr(.x = mda_class,.y = count,.f = function(.x,.y) {
        if(.y < mask_cutoff2) {return(paste0("tiny categories\n(less than ",mask_cutoff2," compounds.)"))} else {return(.x)}
      })
    ) %>%
    group_by(class_less) %>%
    summarise(num = sum(count)) %>%
    setNames(c("mda_class","count"))

  category = list(
    category_minor = category_minor %>% arrange(desc(count)),
    category_major = category_major %>% arrange(desc(count))
  )
  return(category)
}

step10 = percent.tbl %>%
  select(variable_id,Superclass,Class,Subclass) %>%
  distinct() %>%
  setNames(c("compound_id","superclass","class","subclass"))
subclass_catg = mda_sum_categroy2(cate_df = step10,cate_type = "subclass",mask_cutoff = 10,mask_cutoff2 = 3)
dir.create("01.Metbo")

pdf("01.Metbo/03.sub_class_major.pdf",width = 20,height = 15)
pie(subclass_catg$category_major$count,labels = paste0(subclass_catg$category_major$mda_class," (",subclass_catg$category_major$count,")"))
dev.off()
pdf("01.Metbo/03.sub_class_minor.pdf",width = 20,height = 15)
pie(subclass_catg$category_minor$count,labels = paste0(subclass_catg$category_minor$mda_class," (",subclass_catg$category_minor$count,")"))
dev.off()

png("01.Metbo/03.sub_class_major.png",width = 2000,height = 1500,res = 300)
pie(subclass_catg$category_major$count,labels = paste0(subclass_catg$category_major$mda_class," (",subclass_catg$category_major$count,")"))
dev.off()
png("01.Metbo/03.sub_class_minor.png",width = 2000,height = 1500,res = 300)
pie(subclass_catg$category_minor$count,labels = paste0(subclass_catg$category_minor$mda_class," (",subclass_catg$category_minor$count,")"))
dev.off()

sum(subclass_catg$category_minor$count)
sum(subclass_catg$category_major$count)

class_catg = mda_sum_categroy2(cate_df = step10,cate_type = "class",mask_cutoff = 10,mask_cutoff2 = 3)
pdf("01.Metbo/03.class_major.pdf",width = 20,height = 15)
pie(class_catg$category_major$count,labels = paste0(class_catg$category_major$mda_class," (",class_catg$category_major$count,")"))
dev.off()
pdf("01.Metbo/03.class_minor.pdf",width = 20,height = 15)
pie(class_catg$category_minor$count,labels = paste0(class_catg$category_minor$mda_class," (",class_catg$category_minor$count,")"))
dev.off()

png("01.Metbo/03.class_major.png",width = 2000,height = 1500,res = 300)
pie(class_catg$category_major$count,labels = paste0(class_catg$category_major$mda_class," (",class_catg$category_major$count,")"))
dev.off()
png("01.Metbo/03.class_minor.png",width = 2000,height = 1500,res = 300)
pie(class_catg$category_minor$count,labels = paste0(class_catg$category_minor$mda_class," (",class_catg$category_minor$count,")"))
dev.off()
major_tag = subclass_catg$category_major %>%
  mutate(count = 1)

major_tag_class = class_catg$category_major %>%
  mutate(count = 1)
sum(class_catg$category_minor$count)
sum(class_catg$category_major$count)
step10 %>% left_join(major_tag,by = c("subclass" = "mda_class"))%>% select(compound_id,subclass,count) %>% filter(!subclass == "NA") %>%  arrange(count) %>% distinct() -> sub_class_corr_tag
step10 %>% left_join(major_tag_class,by = c("class" = "mda_class")) %>% select(compound_id,class,count) %>%  arrange(count) %>% distinct()-> class_corr_tag
write.csv(step10,"01.Metbo/category.table.csv")
