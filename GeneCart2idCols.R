# tbl has to be a data.frame with IDs as chr and values as numeric
# input table of loadings and 1 and only 1 id col (must be symbol or ensgID)
# only good for mouse and human now - need to add correct getMatch output col names for other species b4 it will work for other species:
GeneCart2idCols=function(tbl=NA,IDcol="",IDtype="",IDspec="",removeDUPgenesymbol=FALSE,removeDUPensembl=TRUE,useNewestVersion=FALSE){

if(IDspec!="human"&IDspec!="mouse"){print("this function only works for mouse or human currently - see function comments for more info.");return(NULL)}

library(SJD)
IDcolNUM=which(colnames(tbl)==IDcol)

if(IDtype=="symbol"){
moreIDs=getMatch(genes=tbl[,IDcolNUM],inSpecies=IDspec,inType=IDtype,newSpecies=IDspec,useNewestVersion=useNewestVersion)
if(IDspec=="human"){tbl2=cbind(moreIDs$Gene.stable.ID.human,tbl[,IDcolNUM],tbl[,-IDcolNUM])}
if(IDspec=="mouse"){tbl2=cbind(moreIDs$Gene.stable.ID.mouse,tbl[,IDcolNUM],tbl[,-IDcolNUM])}
}

if(IDtype=="ensembl"){
moreIDs=getMatch(genes=tbl[,IDcolNUM],inSpecies=IDspec,inType=IDtype,newSpecies=IDspec,useNewestVersion=useNewestVersion)
if(IDspec=="human"){tbl2=cbind(tbl[,IDcolNUM],moreIDs$HGNC.symbol.human,tbl[,-IDcolNUM])}
if(IDspec=="mouse"){tbl2=cbind(tbl[,IDcolNUM],moreIDs$MGI.symbol.mouse,tbl[,-IDcolNUM])}
}

colnames(tbl2)[1:2]=c("ensemblGeneID","GeneSymbol")

# converts dups and NAs and ""s into dumby IDs to prevent loss of rows for which we only have 1 type of ID
#
# ensembl
print("converting dups and NAs in ensemblIDs into dumby ensemblIDs to prevent loss of rows for which we only have gene symbols (for NeMO only):")
indxDROP=duplicated(tbl2[,1])|is.na(tbl2[,1])|tbl2[,1]==""
print(paste0("Ensembl IDs: ",sum(indxDROP)," dups/NAs/empties will be converted into dumby IDs (of ",length(indxDROP)," total Ensembl IDs)."))
seqDUBMYids=paste("NOensemblIDmapped",c(1:sum(indxDROP)),sep=".")
tbl2[,1][indxDROP]=seqDUBMYids
# symbol
print("converting dups and NAs in symbols into dumby gene symbols to prevent loss of rows for which we only have ensembl IDs (for NeMO only):")
indxDROP=duplicated(tbl2[,2])|is.na(tbl2[,2])|tbl2[,2]==""
print(paste0("Gene Symbols: ",sum(indxDROP)," dups/NAs/empties will be converted into dumby IDs (of ",length(indxDROP)," total gene symbols)."))
seqDUBMYids=paste("NOsymbolMapped",c(1:sum(indxDROP)),sep=".")
tbl2[,2][indxDROP]=seqDUBMYids

# print(str(tbl2))
# # print(str(tbl2$GeneSymbol))
# print(str(tbl2[,"GeneSymbol"]))

# ensure genesymbols have no more than 20 characters:
# tbl2$GeneSymbol=sapply(tbl2$GeneSymbol,substr,1,20)
tbl2[,"GeneSymbol"]=sapply(tbl2[,"GeneSymbol"],substr,1,20)

# colnames(tbl2)=gsub(" ","_",colnames(tbl2))
# colnames(tbl2)=gsub("-","_",colnames(tbl2))

if(removeDUPgenesymbol){
tbl2=tbl2[!duplicated(tbl2[,"GeneSymbol"]),]
print(paste0("removeDUPgenesymbol=TRUE; dropped ",sum(duplicated(tbl2[,"GeneSymbol"]))," duplicated GeneSymbols."))
}

if(removeDUPensembl){
tbl2=tbl2[!duplicated(tbl2[,"ensemblGeneID"]),]
print(paste0("removeDUPensembl=TRUE; dropped ",sum(duplicated(tbl2[,"ensemblGeneID"]))," duplicated ensemblGeneIDs."))
}

return(tbl2)
}
