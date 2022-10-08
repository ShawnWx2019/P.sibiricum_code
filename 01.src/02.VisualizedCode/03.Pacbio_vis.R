###########################################################
#      Prj: Multi-omics data analysis of P.sibiricum.
#      Assignment: Upset of annotation
#      Date: Mar 21, 2022
#      Author: Shawn Wang <shawnwang2016@126.com>
#      Location: HENU, Kaifeng, Henan, China
###########################################################
library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
# dir.create("../03.Progress/02.Pacbio/Annotation/",showWarnings = F,recursive = T)
# setwd("../03.Progress/02.Pacbio/Annotation/")
Pacbio_anno <- readxl::read_xlsx("Annotation.xlsx") %>% distinct()
# make a list
Pacbio_anno <-
  Pacbio_anno %>% mutate(across(
   .cols = !GeneID,
   ~ case_when(
     . == "--" ~ 0,
     TRUE ~ 1
   )
  )) %>%
  mutate(sum = apply(.[,-1],1,sum)) %>%
  group_by(GeneID) %>%
  top_n(n = 1,wt = sum) %>%
  select(-sum) %>%
  ungroup()
colnames(Pacbio_anno) = c("Gene_ID","GO","KEGG","KOG","NR","NT","Swissprot")

data_set <- map(.x = colnames(Pacbio_anno)[-1],.f = function(.x) {
  set <-
    Pacbio_anno %>%
    select(Gene_ID,.x) %>%
    dplyr::filter(.[,2] == 1)%>%
    distinct() %>%
    pull(Gene_ID)
})

names(data_set) = colnames(Pacbio_anno)[-1]


data_mat <-
  Pacbio_anno %>%
  column_to_rownames("Gene_ID")

m1 = make_comb_mat(data_set)
m <- m1
ss = set_size(m)
cs = comb_size(m)
comb_degree(m)
ht = UpSet(m,
           set_order = order(ss),comb_col = "#288CD2",
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Annotation Intersections" = anno_barplot(cs,
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE,
                                                  gp = gpar(fill = "#288CD2"),
                                                  height = unit(8,"cm")
             ),
             annotation_name_side = "left",
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Annotation in each database" = anno_barplot(-ss,
                                               baseline = 0,
                                               border = FALSE,
                                               gp = gpar(fill = "#3DFF92"),
                                               width = unit(3, "cm"),
                                               height = unit(3,"cm")
             ),
             set_name = anno_text(set_name(m),
                                  location = 0.5,
                                  just = "center"
                                 # width = max_text_width(set_name(m)) + unit(4, "mm")
                                 )
           ),
           right_annotation = NULL,
           show_row_names = FALSE)
pdf("UpSetPlot.pdf",width = 10,height = 6)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Annotation Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
            default.units = "native", just = c("left", "bottom"),
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
dev.off()

data_mat %>%
  filter(rowSums(.) >2) %>% nrow()

Annotation_all <- read.delim("../../../02.Data/02.Pacbio_raw/Integrated_Function.annotation.xls")
species_percentage <-
Annotation_all %>%
  select(Annotation) %>%
  mutate(Annotation =
           str_extract(Annotation,"(?<=\\[).*(?=\\])"),
         n = 1) %>%
  mutate(Annotation = case_when(
    is.na(Annotation) ~ "Not annotated",
    TRUE ~ Annotation
  )) %>%
  group_by(Annotation) %>%
  summarise(count = sum(n)) %>%
  mutate(percentage = count/143998) %>%
  arrange(desc(count))
writexl::write_xlsx(species_percentage,"../../../04.Result/01.Tables/TableS7.Annotation_Summary.xlsx")

writexl::write_xlsx(Annotation_all,"../../../04.Result/01.Tables/TableS6.Annotation_all.xlsx")

Polygonatum.tbl <-
Annotation_all %>%
  select(Annotation,Identity) %>%
  mutate(Annotation =
           str_extract(Annotation,"(?<=\\[).*(?=\\])"),
         Identity =
           str_extract(Identity,"(?<=\\().*(?=\\))"),
         n = 1) %>%
  mutate(Annotation = case_when(
    is.na(Annotation) ~ "Not annotated",
    TRUE ~ Annotation
  )) %>%
  filter(
    str_detect(Annotation,"Polygonatum")
  ) %>%
  mutate(Identity = as.numeric(Identity)) %>%
  arrange(Identity) %>%
  group_by(Annotation) %>%
  summarise(mean = mean(Identity),
            sum = sum(n)) %>%
  arrange(mean,sum)
colnames(Polygonatum.tbl) = c("Species","Identity","hit number")
writexl::write_xlsx(Polygonatum.tbl,"../../../04.Result/01.Tables/TableS8.xlsx")
Annotation_all %>%
  select(Annotation,Identity) %>%
  mutate(Annotation =
           str_extract(Annotation,"(?<=\\[).*(?=\\])"),
         Identity =
           str_extract(Identity,"(?<=\\().*(?=\\))"),
         n = 1) %>%
  mutate(Annotation = case_when(
    is.na(Annotation) ~ "Not annotated",
    TRUE ~ Annotation
  )) %>%
  filter(
    str_detect(Annotation,"Polygonatum")
  ) %>%
  mutate(Identity = as.numeric(Identity)) %>%
  arrange(Identity) %>%
  pull(Identity) %>% fivenum()

tmp_a <- read.delim("../../../02.Data/02.Pacbio_raw/UniIso.fa.KOG_class.xls")
tmp_a %>%
  select(Functional_categories) %>%
  mutate(n = 1) %>%
  group_by(Functional_categories) %>%
  summarise(sums = sum(n))
tmp_b<- list()
for (i in 1:length(tmp_a$Functional_categories %>% unique())) {
  tmp_b[[i]] = data.frame(
    category = unique(tmp_a$Functional_categories)[i],
    count = tmp_a %>%
      filter(Functional_categories == unique(Functional_categories)[i]) %>%
      select(X.GeneID,Functional_categories) %>%
      pull(X.GeneID) %>% unique() %>% length()
  )
}
tmp_b = bind_rows(tmp_b)
writexl::write_xlsx(tmp_b,"KOG_summary.xlsx")

GO_count <- read.delim("../../../02.Data/02.Pacbio_raw/GO.anno.xls")

head(GO_count)
GO_clean <-
GO_count %>%
  group_by(X.GeneID) %>%
  select(-GOTerms_num) %>%
  pivot_longer(!X.GeneID,names_to = "tag",values_to = "term") %>%
  filter(term != "") %>%
  select(-tag) %>%
  ungroup() %>%
  mutate(
   term = gsub("\\(biological_process","|biological_process",term) %>%
     gsub("\\(molecular_function","|molecular_function",term) %>%
     gsub("\\(cellular_component","|cellular_component",term)
  ) %>%
  mutate(
    GO_ID = str_split(term,"\\: ",n = 2,simplify = T)[,1],
    Description = str_extract(term,"(?<=\\: ).*(?= \\()"),
    category = str_extract(term,"(?<=\\|).*(?=\\))")
  ) %>%
  select(-term) %>%
  distinct()

T2G = GO_clean %>%
  select(GO_ID,X.GeneID) %>%
  setNames(c("TERM","GENE"))
T2N = GO_clean %>%
  select(GO_ID,Description) %>%
  setNames(c("TERM","NAME"))
GO_clean %>%
  select(X.GeneID,category) %>%
  group_by(category) %>%
  mutate(n = 1) %>%
  summarise(count = sum(n))


# kegg --------------------------------------------------------------------

kegg_category <- read.delim("/Users/shawn/SynologyDrive/Project/11.HJ/01.Data/report/result/Annotation/HJ/KEGG/KEGG_classification_count.txt")
kegg_category %>%
  group_by(X.Pathway.Hierarchy1) %>%
  summarise(n = sum(Gene.Number))


# CDS length --------------------------------------------------------------
cd_hit <- read.delim("/Users/shawn/SynologyDrive/Project/11.HJ/01.Data/report/result/CD-hit-est/HJ/cd_hit_id.xls",header = F)
cds_info <- read.delim("/Users/shawn/SynologyDrive/Project/11.HJ/01.Data/report/result/Transdecoder/HJ/transdecoder.list.xls",header = F)

cds.length <- read.delim("/Users/shawn/SynologyDrive/Project/11.HJ/01.Data/report/result/Transdecoder/HJ/length.xls",header = F)

cds_length_histogram <-
cds.length %>%
  select(V2) %>%
  ggplot(aes(V2)) + geom_histogram(
    stat = 'bin',fill = "blue",color = 'black',alpha = .6,binwidth = 300
  )+
  xlab("Length (nt)")+
  ylab("")+
  ggtitle("CDS sequence length distribution")+
  theme_bw()+
  theme(
    rect = element_rect(size = 1.5,colour = 'black'),
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title = element_text(size = 13,colour = 'black'),
    title = element_text(size = 13,colour = 'black',face = 'bold')
  )

lncRNA_id <- read.delim("/Users/shawn/SynologyDrive/Project/11.HJ/01.Data/report/result/LncRNA/HJ/HJ_lncrna_id.txt",header = F)

cds.id <- cds_info$V1 %>% str_extract("(?<=\\>)HJ_1-10k_transcript/\\d++") %>% unique()

lnc.id <- lncRNA_id$V1 %>% str_extract("(?<=\\>).*(?= f)") %>% unique()

cd_hit.id = cd_hit$V1 %>% str_extract("(?<=\\>)HJ_1-10k_transcript/\\d++") %>% unique()
v_list = list(
  CDS = cds.id,
  LncRNA = lnc.id,
  `non-redundant` = cd_hit.id
)

venn_plt <-
ggvenn(v_list)

lnc.length <- read.delim("/Users/shawn/SynologyDrive/Project/11.HJ/01.Data/report/result/LncRNA/HJ/lnc_length.xls",header = F)

lnc_length_histogram <-
  lnc.length %>%
  select(V2) %>%
  ggplot(aes(V2)) + geom_histogram(
    stat = 'bin',fill = "blue",color = 'black',alpha = .6,binwidth = 300
  )+
  xlab("Length (nt)")+
  ylab("")+
  ggtitle("lncRNA sequence length distribution")+
  theme_bw()+
  theme(
    rect = element_rect(size = 1.5,colour = 'black'),
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title = element_text(size = 14,colour = 'black'),
    title = element_text(size = 13,colour = 'black',face = 'bold')
  )

figureS4 <-cds_length_histogram + venn_plt + lnc_length_histogram + plot_annotation(tag_levels = "A")

ggsave(filename = "../../../04.Result/02.Figures/FigS4.non_redundant.pdf",
       plot = figureS4,
       width = 24,
       height = 8)

ggsave(filename = "../../../05.Report/FigS4.non_redundant.png",
       plot = figureS4,
       width = 24,
       height = 8)

