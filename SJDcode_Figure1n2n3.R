#####################################################################################################################################################################
# # wgets from nemo for fig 1/2/3: from unix command line
#
cd /path/
# lamanno:
wget -O LaManno.h5ad 'https://nemoanalytics.org/cgi/download_source_file.cgi?type=h5ad&share_id=dfda4722'
# micali: 
wget -O Micali.h5ad 'https://nemoanalytics.org/cgi/download_source_file.cgi?type=h5ad&share_id=e9947c2c'
# trevino: 
wget -O Trevino.h5ad 'https://nemoanalytics.org/cgi/download_source_file.cgi?type=h5ad&share_id=b63f99b1'

##############################################################################################################################################################
# # begin R
library(anndata)
musmus <- read_h5ad("/path/LaManno.h5ad")# La Manno 2021
macmul <- read_h5ad("/path/Micali.h5ad")# Micali 2023
homsap <- read_h5ad("/path/Trevino.h5ad")# Trevino 2021

###############################################################################
# assemble data for joint decomposition
# library(devtools)
# install_github("CHuanSite/SJD")
library(SJD)

baseDIR="/path/"
lbb1="MidGestExcitNroGen.scRNAseq.Mm1MacMul1Hs1.TEST"
META.list=list()
EXPRS.list=list()

META.list$laManno=cbind(colnames(t(as.matrix(musmus$X))),musmus$obs)
META.list$Micali=cbind(colnames(t(as.matrix(macmul$X))),macmul$obs)
META.list$Trevino=cbind(colnames(t(as.matrix(homsap$X))),homsap$obs)

EXPRS.list$laManno=t(as.matrix(musmus$X))
EXPRS.list$Micali=t(as.matrix(macmul$X))
EXPRS.list$Trevino=t(as.matrix(homsap$X))

rm(musmus)
rm(macmul)
rm(homsap)
gc()

# str(META.list)
names(META.list)
# [1] "laManno" "Micali"  "Trevino"

str(EXPRS.list)
# List of 3
 # $ laManno: num [1:31053, 1:38668] 14.05 8.57 11.02 10.15 9.15 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:31053] "ENSMUSG00000092341" "ENSMUSG00000049775" "ENSMUSG00000050708" "ENSMUSG00000001525" ...
  # .. .. ..- attr(*, "name")= chr "gene"
  # .. ..$ : chr [1:38668] "X10X70_5_A_1.AGGAATGATTTCACx" "X10X32_2_A_1.TTATGAGAGATGAAx" "X10X13_4_A_1.CATTACACCTAAGCx" "X10X73_3_A_1.AAGGTCTGCCAATGx" ...
  # .. .. ..- attr(*, "name")= chr "observations"
 # $ Micali : num [1:34619, 1:70407] 0 0 0 0 0 0 0 0 0 0 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:34619] "ENSG00000185220" "NOensemblIDmapped.1" "NOensemblIDmapped.2" "ENSG00000171163" ...
  # .. .. ..- attr(*, "name")= chr "gene"
  # .. ..$ : chr [1:70407] "E54_FR_1231_AAACCTGAGAACTGTA" "E54_FR_1231_AAACCTGAGCCATCGC" "E54_FR_1231_AAACCTGAGCGAGAAA" "E54_FR_1231_AAACCTGAGGATGGAA" ...
  # .. .. ..- attr(*, "name")= chr "observations"
 # $ Trevino: num [1:33325, 1:43340] 0 0 0 0 0 0 0 0 0 0 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:33325] "NOensemblIDmapped.1" "ENSG00000237613" "ENSG00000186092" "NOensemblIDmapped.2" ...
  # .. .. ..- attr(*, "name")= chr "gene"
  # .. ..$ : chr [1:43340] "hft_w20_p3_r1_AAACCCACATAGTCAC" "hft_w20_p3_r1_AAACCCAGTACAGGTG" "hft_w20_p3_r1_AAACCCAGTACGGTTT" "hft_w20_p3_r1_AAACCCAGTACTCGCG" ...
  # .. .. ..- attr(*, "name")= chr "observations"

############################
# info about META.sample for automated plotting post-SJD
SampleMetaNamesTable=data.frame(
row.names = names(EXPRS.list),
Type=c("2Dscatter","2Dscatter","2Dscatter"),#'Yaxis'
XaxisColumn=c("UMAP_1","UMAP_1","UMAP_0"),
YaxisColumn=c("UMAP_2","UMAP_2","UMAP_1"),
COLaxisColumn=c("SJDscores","SJDscores","SJDscores"),
PCHColumn=c("","",""),
cexx=c(1,1,1)
)
print(SampleMetaNamesTable)
             # Type XaxisColumn YaxisColumn COLaxisColumn PCHColumn cexx
# laManno 2Dscatter      UMAP_1      UMAP_2     SJDscores              1
# Micali  2Dscatter      UMAP_1      UMAP_2     SJDscores              1
# Trevino 2Dscatter      UMAP_0      UMAP_1     SJDscores              1

###############################
# orthologue mapping and data prep for SJD

SJDdataIN=sjdWrap(
data.list=EXPRS.list,
species.vector=c("mouse","human","human"),
geneType.vector=c("ensembl","ensembl","ensembl"),
geneType.out="ensembl",species.out="human"
)
# Using biomaRt to connect gene IDs across 3 datasets:
# Getting biomaRt IDs for dataset 1
# You have input  31053  genes
# We found  23449  matches
# 4609  of those are duplicates and only keeping the 1st of each
# Getting biomaRt IDs for dataset 2
# You have input  34619  genes
# We found  17729  matches
# 0  of those are duplicates and only keeping the 1st of each
# Getting biomaRt IDs for dataset 3
# You have input  33325  genes
# We found  21518  matches
# 0  of those are duplicates and only keeping the 1st of each
# constructed 3 tables of cross-species matching genes
# we found 15625 shared genes in 3 datasets
# new data list of 3 datasets constructed

str(SJDdataIN)
lapply(SJDdataIN,dim)
# $laManno
# [1] 15625 38668
# $Micali
# [1] 15625 70407
# $Trevino
# [1] 15625 43340

rm(EXPRS.list)

####################################################################################
# searhcing for common elements of variation across all 3 matrices
grp=list(Shared.Mammal=c(1:3))

########
# 7 dimension decomposition
dims=c(7)# must have same length as "grp"
lbb2=paste(lbb1,"p",dims[1],sep=".")
print(Sys.time())
jointNMF=jointNMF(dataset=SJDdataIN,group=grp,comp_num=dims)
print(Sys.time())
save(SampleMetaNamesTable,grp,dims,META.list,jointNMF,file=paste(baseDIR,"SJDdataOUT_",lbb2,"_jointNMF_fullShareONLY.RData",sep=""))
source("/path/GeneCart2idCols.R")
tbl2=GeneCart2idCols(tbl=data.frame(ensembl=rownames(jointNMF$linked_component_list$Shared.Mammal),jointNMF$linked_component_list$Shared.Mammal),IDcol="ensembl",IDtype="ensembl",IDspec="human",useNewestVersion=TRUE)
str(tbl2)
write.table(tbl2,row.names=FALSE,col.names=TRUE,sep="\t",file=paste0(baseDIR,lbb2,"_genelist.tab"))# this file can be uploaded to NeMO Analytics to create a gene list for projection
################
# SJD plotting
library(ggplot2)
library(gridExtra)
alg="jointNMF"
scrs=jointNMF$score_list;grpp=names(scrs[[1]])[1];kt=dim(scrs[[1]][[grpp]])[1];kk=1:kt;Nrows=length(kk)
Ncols=length(scrs)
SJDScorePlotter.obj=SJDScorePlotter(SJDalg=alg,scores=scrs,lbb=lbb2,info=META.list,SampleMetaNamesTable=SampleMetaNamesTable)
assemble.byComponent=assemble.byComponent(SJDScorePlotter.obj=SJDScorePlotter.obj,component=kk,SJD_algorithm=alg,group='Shared.Mammal')
adj=5
pdf(width=Ncols*adj,height=Nrows*adj,paste0(baseDIR,"SJDdataOUT_",lbb2,"_",alg,"_",grpp,"_k",kk[1],".to.",kk[length(kk)],"of",kt,"new.pdf"))
g = grid.arrange(grobs=assemble.byComponent,nrow=Nrows,ncol=Ncols)
dev.off()


########
# 40 dimension decomposition
dims=c(40)# must have same length as "grp"
lbb2=paste(lbb1,"p",dims[1],sep=".")
print(Sys.time())
jointNMF=jointNMF(dataset=SJDdataIN,group=grp,comp_num=dims)
print(Sys.time())
save(SampleMetaNamesTable,grp,dims,META.list,jointNMF,file=paste(baseDIR,"SJDdataOUT_",lbb2,"_jointNMF_fullShareONLY.RData",sep=""))
source("/path/GeneCart2idCols.R")
tbl2=GeneCart2idCols(tbl=data.frame(ensembl=rownames(jointNMF$linked_component_list$Shared.Mammal),jointNMF$linked_component_list$Shared.Mammal),IDcol="ensembl",IDtype="ensembl",IDspec="human",useNewestVersion=TRUE)
str(tbl2)
write.table(tbl2,row.names=FALSE,col.names=TRUE,sep="\t",file=paste0(baseDIR,lbb2,"_genelist.tab"))# this file can be uploaded to NeMO Analytics to create a gene list for projection
################
# SJD plotting
library(ggplot2)
library(gridExtra)
alg="jointNMF"
scrs=jointNMF$score_list;grpp=names(scrs[[1]])[1];kt=dim(scrs[[1]][[grpp]])[1];kk=1:kt;Nrows=length(kk)
Ncols=length(scrs)
SJDScorePlotter.obj=SJDScorePlotter(SJDalg=alg,scores=scrs,lbb=lbb2,info=META.list,SampleMetaNamesTable=SampleMetaNamesTable)
assemble.byComponent=assemble.byComponent(SJDScorePlotter.obj=SJDScorePlotter.obj,component=kk,SJD_algorithm=alg,group=grpp)
adj=5
pdf(width=Ncols*adj,height=Nrows*adj,paste0(baseDIR,"SJDdataOUT_",lbb2,"_",alg,"_",grpp,"_k",kk[1],".to.",kk[length(kk)],"of",kt,"new.pdf"))
g = grid.arrange(grobs=assemble.byComponent,nrow=Nrows,ncol=Ncols)
dev.off()



#####################
quit(save='no')


