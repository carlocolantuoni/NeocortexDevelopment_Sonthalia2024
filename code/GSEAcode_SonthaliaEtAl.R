library(SJD)
library(fgsea)
library(limma)
library(dplyr)

comps=read.delim(as.is=TRUE,file="SonthaliaEtAl_SupplTable1_p7CtxDev_GeneLoadings.txt")
rownames(comps)=comps$GeneSymbol
comps=comps[,3:9]
str(comps)

# load gene groups -  files from: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/
go<-  gmtPathways("/path/c5.go.v2023.1.Hs.symbols.gmt.txt")
hallmark <-  gmtPathways("/path/h.all.v2023.1.Hs.symbols.gmt.txt")
kegg <-  gmtPathways("/path/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt")
tft <-  gmtPathways("/path/c3.tft.v2023.1.Hs.symbols.gmt.txt")
load("/path/CtxDiseaseGeneListsPLUS.RData")# from https://pubmed.ncbi.nlm.nih.gov/40634286/
names(CtxDiseaseGeneLists)=paste('DiseaseLists',names(CtxDiseaseGeneLists),sep='_')
names(tft)=paste('TFtargets',names(tft),sep='_')

all_sets=c(go,hallmark,kegg,tft,CtxDiseaseGeneLists)

sizes=c(lengths(go),lengths(hallmark),lengths(kegg),lengths(tft),lengths(CtxDiseaseGeneLists))

intersection_sizes=c()
for (i in all_sets){
  intersection_sizes=c(intersection_sizes,(sum(rownames(comps) %in% i)))
}

geneset_db=c(rep('GO',times=length(go)),rep('Hallmark',times=length(hallmark)),rep('Kegg',times=length(kegg)),rep('TFT',times=length(tft)),rep('DiseaseLists',times=length(CtxDiseaseGeneLists)))
table(geneset_db)
# DiseaseLists           GO     Hallmark         Kegg          TFT
         # 106        10532           50          186         1115

gst=function(subcomp,pattern,sets){
  pvals_subcomp=c()
  for (i in sets){pvals_subcomp=c(pvals_subcomp, wilcoxGST(which(rownames(pattern) %in% i), pattern[,subcomp],type='f',alternative='mixed'))}
  return(pvals_subcomp)}

rng=(1:dim(comps)[2])

subcomp_ALL=data.frame(sapply(rng,gst,pattern=comps,sets=all_sets))

subcomp_ALL=cbind(geneset_db,sizes,intersection_sizes,names(all_sets),subcomp_ALL)
colnames(subcomp_ALL)=c("GeneSetDB","size","intersection","GeneSet",paste0("pattern",1:dim(comps)[2]))

write.csv(subcomp_ALL,file = "SonthaliaEtAl_SupplTable1_p7CtxDev_GST.csv")

