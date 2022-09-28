# PID of current job: 2998577
mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c("C21589","C21512","C21136","C21044","C18402","C17582","C17392","C17222","C16357","C12254","C11644","C11073","C10897","C10435","C07559","C07207","C06435","C06427","C06269","C06218","C06217","C06082","C05662","C05488","C04742","C03852","C03274","C02217","C02140","C02132","C01953","C01620","C01595","C01324","C01187","C01157","C00864","C00860","C00624","C00492","C00811","C07189","C12204","C04738","C01026","C00383","C00160","C00094","C01594","C00299","C08328","C01709","C09638","C06470","C00601","C10320","C08729","C08885")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "kegg");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetKEGG.PathLib(mSet, "ath", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
mSet<-PlotPathSummary(mSet, F, "path_view_1_", "png", 72, width=NA, NA, NA )
mSet<-PlotKEGGPath(mSet, "Biosynthesis of unsaturated fatty acids",576, 480, "png", NULL)
mSet<-RerenderMetPAGraph(mSet, "zoom1664116384447.png",576.0, 480.0, 100.0)
mSet<-PlotKEGGPath(mSet, "Linoleic acid metabolism",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Biosynthesis of secondary metabolites - unclassified",576, 480, "png", NULL)
