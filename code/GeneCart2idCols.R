# tbl has to be a data.frame with IDs as chr and values as numeric
# input table of loadings and 1 and only 1 id col (must be symbol or ensgID)
#
# eg below : 
# tbl0 is a data.frame with column "symbols" containing gene symbols and the other columns containing numeric weights
# tbl2 is the outpt written to a file that can be uploaded to nemo analytics
#
# source(/path/GeneCart2idCols.R)
# tbl2=GeneCart2idCols(tbl=tbl0,IDcol="symbols",IDtype="symbol",IDspec="human",useNewestVersion=TRUE)
# str(tbl2)
# write.table(tbl2,row.names=FALSE,col.names=TRUE,sep="\t",file=paste0("/path/file.tab"))

GeneCart2idCols=function(tbl=NA,IDcol="gene",IDtype="symbol",IDspec="human",removeDUPgenesymbol=FALSE,removeDUPensembl=TRUE,useNewestVersion=TRUE){

library(SJD)
IDcolNUM=which(colnames(tbl)==IDcol)

if(IDtype=="symbol"){
moreIDs=getMatch(genes=tbl[,IDcolNUM],inSpecies=IDspec,inType=IDtype,newSpecies=IDspec,useNewestVersion=useNewestVersion)
tbl2=cbind(moreIDs$ensembl_gene_id,tbl[,IDcolNUM],tbl[,-IDcolNUM])
}

if(IDtype=="ensembl"){
moreIDs=getMatch(genes=tbl[,IDcolNUM],inSpecies=IDspec,inType=IDtype,newSpecies=IDspec,useNewestVersion=useNewestVersion)
tbl2=cbind(tbl[,IDcolNUM],moreIDs$external_gene_name,tbl[,-IDcolNUM])
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

