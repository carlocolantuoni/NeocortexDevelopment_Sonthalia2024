#####################################################################################################################################################################
# # wgets from nemo for fig 4/5/6: from unix command line
#
cd /path/
# jorstad 2023: 
wget https://nemoanalytics.org/datasets/e5dcba0a-07e5-c774-2f58-279acb74e1c7.h5ad
# bakken 2021: 
wget https://nemoanalytics.org/datasets/69520ae8-8b1d-2a16-1118-9be6b66d4c58.h5ad

##############################################################################################################################################################
# # begin R
library(anndata)
jor <- read_h5ad("/path/e5dcba0a-07e5-c774-2f58-279acb74e1c7.h5ad")
bak <- read_h5ad("/path/69520ae8-8b1d-2a16-1118-9be6b66d4c58.h5ad")

#indices to later split into individal donor matrices
str(jor$obs)
indxJOR1=jor$obs$external_donor_name_label=='H200.1025';length(indxJOR1);sum(indxJOR1)#31229, 6658
indxJOR2=jor$obs$external_donor_name_label=='H200.1030';length(indxJOR2);sum(indxJOR2)#31229, 12354
indxJOR3=jor$obs$external_donor_name_label=='H200.1023';length(indxJOR3);sum(indxJOR3)#31229, 12217
str(bak$obs)
indxBAK1=bak$obs$donor_id=='H18.30.001';length(indxBAK1);sum(indxBAK1)#7805, 4576
indxBAK2=bak$obs$donor_id=='H18.30.002';length(indxBAK2);sum(indxBAK2)#7805, 3229

###############################################################################
# assemble data for joint decomposition
# library(devtools)
# install_github("CHuanSite/SJD")
library(SJD)

baseDIR="/path/"
lbb1="AdultNeoctxNeuron.snRNAseq.JorBak.TEST"
META.list=list()
EXPRS.list=list()

META.list$jor1=cbind(colnames(t(as.matrix(jor$X[indxJOR1,]))),jor$obs[indxJOR1,])
META.list$jor2=cbind(colnames(t(as.matrix(jor$X[indxJOR2,]))),jor$obs[indxJOR2,])
META.list$jor3=cbind(colnames(t(as.matrix(jor$X[indxJOR3,]))),jor$obs[indxJOR3,])
META.list$bak1=cbind(colnames(t(as.matrix(bak$X[indxBAK1,]))),bak$obs[indxBAK1,])
META.list$bak2=cbind(colnames(t(as.matrix(bak$X[indxBAK2,]))),bak$obs[indxBAK2,])

EXPRS.list$jor1=t(as.matrix(jor$X[indxJOR1,]))
EXPRS.list$jor2=t(as.matrix(jor$X[indxJOR2,]))
EXPRS.list$jor3=t(as.matrix(jor$X[indxJOR3,]))
EXPRS.list$bak1=t(as.matrix(bak$X[indxBAK1,]))
EXPRS.list$bak2=t(as.matrix(bak$X[indxBAK2,]))

rm(jor)
rm(bak)
gc()

# str(META.list)
names(META.list)
# [1] "jor1" "jor2" "jor3" "bak1" "bak2"

str(EXPRS.list)
# List of 5
 # $ jor1: num [1:50281, 1:6658] 0 0 0 0 0 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:50281] "NOensemblIDmapped.1" "NOensemblIDmapped.2" "NOensemblIDmapped.3" "NOensemblIDmapped.4" ...
  # .. ..$ : chr [1:6658] "F2S4_160113_028_B01" "F2S4_160113_028_C01" "F2S4_160113_029_D01" "F2S4_160113_030_A01" ...
 # $ jor2: num [1:50281, 1:12354] 0 0 0 0 0 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:50281] "NOensemblIDmapped.1" "NOensemblIDmapped.2" "NOensemblIDmapped.3" "NOensemblIDmapped.4" ...
  # .. ..$ : chr [1:12354] "F1S4_160106_001_D01" "F1S4_160106_001_E01" "F1S4_160106_001_G01" "F1S4_160106_001_H01" ...
 # $ jor3: num [1:50281, 1:12217] 0 0 0 0 0 0 0 0 0 0 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:50281] "NOensemblIDmapped.1" "NOensemblIDmapped.2" "NOensemblIDmapped.3" "NOensemblIDmapped.4" ...
  # .. ..$ : chr [1:12217] "F1S4_161026_001_B01" "F1S4_161026_001_F01" "F1S4_161026_002_C01" "F1S4_161026_076_A01" ...
 # $ bak1: num [1:14544, 1:4576] 0.94 0 0 0 0 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:14544] "ENSG00000188976" "ENSG00000187583" "ENSG00000187642" "ENSG00000187608" ...
  # .. ..$ : chr [1:4576] "AAACCCAAGTATGGCG.21L8TX_180927_001_A01" "AAACGAAGTATGAAGT.21L8TX_180927_001_A01" "AAACGCTAGTGCACCC.21L8TX_180927_001_A01" "AAAGGATAGAATTCAG.21L8TX_180927_001_A01" ...
 # $ bak2: num [1:14544, 1:3229] 0 0 0 0 0 ...
  # ..- attr(*, "dimnames")=List of 2
  # .. ..$ : chr [1:14544] "ENSG00000188976" "ENSG00000187583" "ENSG00000187642" "ENSG00000187608" ...
  # .. ..$ : chr [1:3229] "AAACCCACAGCTGTAT.25L8TX_180927_003_A01" "AAACCCACATTCAGCA.25L8TX_180927_003_A01" "AAACCCATCTAACGGT.25L8TX_180927_003_A01" "AAACGAAAGGAAGTCC.25L8TX_180927_003_A01" ...

# SampleMetaNamesTable=data.frame(row.names = names(EXPRS.list),
	# Type=c('Yaxis','Yaxis','Yaxis','Yaxis','Yaxis'),#'2Dscatter'
	# XaxisColumn=c("cortical_layer_label","cortical_layer_label","cortical_layer_label","cross_species_cluster_label","cross_species_cluster_label"),
	# YaxisColumn=c("SJDscores","SJDscores","SJDscores","SJDscores","SJDscores"),
	# PCHColumn=c("","","","",""),
	# COLaxisColumn=c("cortical_layer_color","cortical_layer_color","cortical_layer_color","cross_species_cluster_color","cross_species_cluster_color"),
	# cexx=c(1,1,1,1,1))

SampleMetaNamesTable=data.frame(row.names = names(EXPRS.list),
	Type=c('2Dscatter','2Dscatter','2Dscatter','2Dscatter','2Dscatter'),#'Yaxis'
	XaxisColumn=c("tsne_1","tsne_1","tsne_1","tSNE_1","tSNE_1"),
	YaxisColumn=c("tsne_2","tsne_2","tsne_2","tSNE_2","tSNE_2"),
	PCHColumn=c("","","","",""),
	COLaxisColumn=c("SJDscores","SJDscores","SJDscores","SJDscores","SJDscores"),
	cexx=c(1,1,1,1,1))

SJDdataIN=sjdWrap(
  data.list=EXPRS.list,
  species.vector=c("human","human","human","human","human"),
  geneType.vector=c("ensembl","ensembl","ensembl","ensembl","ensembl"),
  geneType.out="ensembl",species.out="human")

str(SJDdataIN)# -  - original decomp had 14272 gene across 3 matrices - this has 14088
lapply(SJDdataIN,dim)
# $jor1
# [1] 14088  6658
# $jor2
# [1] 14088 12354
# $jor3
# [1] 14088 12217
# $bak1
# [1] 14088  4576
# $bak2
# [1] 14088  3229

rm(EXPRS.list)

####################################################################################
# searhcing for common elements of variation across all 5 matrices
grp=list(shared_all=c(1:5))

# 20 dimension decomposition
dims=c(20) # must have same length as "grp"
lbb2=paste(lbb1,".p",dims[1],sep="")

print(Sys.time())
jointNMF=jointNMF(dataset=SJDdataIN,group=grp,comp_num=dims)
print(Sys.time())
save(SampleMetaNamesTable,grp,dims,META.list,jointNMF,file=paste(baseDIR,"SJDdataOUT_",lbb2,"_jointNMF_fullShareONLY.RData",sep=""))
source("/path/GeneCart2idCols.R")
tbl2=GeneCart2idCols(tbl=data.frame(ensembl=rownames(jointNMF$linked_component_list$Shared.Mammal),jointNMF$linked_component_list$Shared.Mammal),IDcol="ensembl",IDtype="ensembl",IDspec="human",useNewestVersion=TRUE)
str(tbl2)
write.table(tbl2,row.names=FALSE,col.names=TRUE,sep="\t",file=paste0(baseDIR,lbb2,"_genelist.tab"))# this file can be uploaded to NeMO Analytics to create a gene list for projection
#################
# SJD plotting
library(ggplot2)
library(gridExtra)
alg="jointNMF"
scrs=jointNMF$score_list;grp=names(scrs[[1]])[1];kt=dim(scrs[[1]][[grp]])[1];kk=1:kt;Nrows=length(kk)
Ncols=length(scrs)
SJDScorePlotter.obj=SJDScorePlotter(SJDalg=alg,scores=scrs,lbb=lbb2,info=META.list,SampleMetaNamesTable=SampleMetaNamesTable)
assemble.byComponent=assemble.byComponent(SJDScorePlotter.obj=SJDScorePlotter.obj,component=kk,SJD_algorithm=alg,group='shared_all')
adj=5
pdf(width=Ncols*adj,height=Nrows*adj,paste0(baseDIR,"SJDdataOUT_",lbb2,"_",alg,"_",grp,"_k",kk[1],".to.",kk[length(kk)],"of",kt,"new.pdf"))
g = grid.arrange(grobs=assemble.byComponent,nrow=Nrows,ncol=Ncols)
dev.off()

#####################
quit(save='no')

