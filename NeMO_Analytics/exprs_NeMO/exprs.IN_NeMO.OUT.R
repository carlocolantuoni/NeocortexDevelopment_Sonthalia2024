exprs.IN_NeMO.OUT=function(
sampLAB,# tag for filenames
exprs,#rows=genes;columns=samples; fully processed expression data, e.g. log2CPM
cellMETA,#rows=samples;columns=info
geneMETA,#rows=genes;columns=info
geneMETA.symbol.column="GeneSymbol",geneMETA.ensembl.column="ensemblGeneID",#colnames in "geneMETA" arg above; need both
baseDIR="/dcs05/carlo/legacy-dcl01-ccolantu/data/Explr/",#output dir
dropDUPna.flag=FALSE,
cnvrtDUPna.flag=TRUE,
nemoMETA.title="",
nemoMETA.summary="",
nemoMETA.dataset_type="scRNA-seq",
nemoMETA.annotation_source="Ensembl",
nemoMETA.annotation_release_number="103",
nemoMETA.geo_accession="",
nemoMETA.contact_email="",
nemoMETA.contact_institute="",
nemoMETA.contact_name="",
nemoMETA.sample_taxid=9606,# 10090 for mouse
nemoMETA.sample_organism="Homo sapiens",# Mus musculus
nemoMETA.platform_id="",
nemoMETA.instrument_model="",
nemoMETA.library_selection="",
nemoMETA.library_source="",
nemoMETA.library_strategy="",
nemoMETA.units="",
nemoMETA.pubmed_id="",
nemoMETA.tags="tissue, cells, experiment",
NeMO.flag=TRUE,TARball.flag=TRUE
)

{

# no spaces in the sampLAB, otherwise we get issues with tar ball cmd below (and possible id and filename issues elsewhere?)
sampLAB=gsub(" ","",sampLAB)

print("************************************************************************************************")
print("************************************************************************************************")
print("Begin.")
print(sampLAB)
print(Sys.time())
print("************************************************************************************************")
print("************************************************************************************************")

# chk if title will be too long for nemo:
if(NeMO.flag&nchar(paste(nemoMETA.title,"[",sampLAB,"]",sep=""))>255){
print(paste0("length of final pasted NeMO title is too long (>255 characters). Ending processing."))
print(paste(nemoMETA.title,"[",sampLAB,"]",sep=""))
return(NULL)
}

#####################
# add on for NeMO
if(NeMO.flag){

print("**********")
print("NeMO output requested.")
print(sampLAB)
print(Sys.time())

print("**********")
print(paste0("dim(exprs):",dim(exprs)))
print(paste0("range(exprs):",range(exprs)))
print(Sys.time())

# format for NeMO

# dump dups and NAs from geneMETA[,geneMETA.ensembl.column] and geneMETA[,geneMETA.symbol.column]
if(dropDUPna.flag){
# geneMETA.ensembl.column
print("dropping dups and NAs from geneMETA.ensembl.column and data! (for NeMO only):")
indxDROP=duplicated(geneMETA[,geneMETA.ensembl.column])|is.na(geneMETA[,geneMETA.ensembl.column])|geneMETA[,geneMETA.ensembl.column]==""
print(paste0("Ensembl IDs: ",sum(indxDROP)," dups/NAs/empties dropped, of ",length(indxDROP)," total Ensembl IDs."))
exprs=exprs[!indxDROP,]
geneMETA=geneMETA[!indxDROP,]
# geneMETA.symbol.column
print("dropping dups and NAs from geneMETA.symbol.column and data! (for NeMO only):")
indxDROP=duplicated(geneMETA[,geneMETA.symbol.column])|is.na(geneMETA[,geneMETA.symbol.column])|geneMETA[,geneMETA.symbol.column]==""
print(paste0("Gene Symbols: ",sum(indxDROP)," dups/NAs/empties dropped, of ",length(indxDROP)," total gene symbols."))
exprs=exprs[!indxDROP,]
geneMETA=geneMETA[!indxDROP,]
}

# converts dups and NAs and ""s from geneMETA[,geneMETA.ensembl.column] and geneMETA[,geneMETA.symbol.column] into dumby IDs to prevent loss of rows for which we only have 1 type of ID
if(cnvrtDUPna.flag){
if(dropDUPna.flag){print("ERROR: can not request both dropDUPna AND cnvrtDUPna");return()}
# geneMETA.ensembl.column
print("converting dups and NAs from geneMETA.ensembl.column into dumby ensembl gene IDs to prevent loss of rows for which we only have gene symbols (for NeMO only):")
indxDROP=duplicated(geneMETA[,geneMETA.ensembl.column])|is.na(geneMETA[,geneMETA.ensembl.column])|geneMETA[,geneMETA.ensembl.column]==""
print(paste0("Ensembl IDs: ",sum(indxDROP)," dups/NAs/empties will be converted into dumby IDs (of ",length(indxDROP)," total Ensembl IDs)."))
seqDUBMYids=paste("NOensemblIDmapped",c(1:sum(indxDROP)),sep=".")
geneMETA[,geneMETA.ensembl.column][indxDROP]=seqDUBMYids
# geneMETA.symbol.column
print("converting dups and NAs from geneMETA.symbol.column into dumby gene symbols to prevent loss of rows for which we only have ensembl IDs (for NeMO only):")
indxDROP=duplicated(geneMETA[,geneMETA.symbol.column])|is.na(geneMETA[,geneMETA.symbol.column])|geneMETA[,geneMETA.symbol.column]==""
print(paste0("Gene Symbols: ",sum(indxDROP)," dups/NAs/empties will be converted into dumby IDs (of ",length(indxDROP)," total gene symbols)."))
seqDUBMYids=paste("NOsymbolMapped",c(1:sum(indxDROP)),sep=".")
geneMETA[,geneMETA.symbol.column][indxDROP]=seqDUBMYids
}

print("**********")
print(paste0("dim(exprs):",dim(exprs)))
print(paste0("range(exprs):",range(exprs)))
print(Sys.time())

##############
# DATA MTX
print("**********")
print("exprs file")
print(Sys.time())
rownames(exprs)=geneMETA[,geneMETA.ensembl.column]
dataRowNames=rownames(exprs)
exprs=data.frame(cbind(dataRowNames,exprs),stringsAsFactors=FALSE)
colnames(exprs)[1]=""
dataColNames=colnames(exprs)[-1]
print("exprs file out")
print(paste(baseDIR,sampLAB,"_DataMTX.tab",sep=""))
write.table(exprs,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE,file=paste(baseDIR,sampLAB,"_DataMTX.tab",sep=""))

##############
# ROW META
print("**********")
print("geneMETA file")
print(Sys.time())
geneMETA=data.frame(cbind(dataRowNames,geneMETA[,geneMETA.symbol.column]),stringsAsFactors=FALSE)
colnames(geneMETA)=c("gene","gene_symbol")
print("geneMETA file out")
print(paste(baseDIR,sampLAB,"_ROWmeta.tab",sep=""))
write.table(geneMETA,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE,file=paste(baseDIR,sampLAB,"_ROWmeta.tab",sep=""))

print(paste0("dim(geneMETA):",dim(geneMETA)))
print(Sys.time())

##############
# CELL/COL META
print("**********")
print("cellMETA file")
print(Sys.time())
cellMETA=data.frame(cbind(dataColNames,cellMETA),stringsAsFactors=FALSE)
colnames(cellMETA)[1]="observations"
print("cellMETA file out")
print(paste(baseDIR,sampLAB,"_COLmeta.tab",sep=""))
write.table(cellMETA,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE,file=paste(baseDIR,sampLAB,"_COLmeta.tab",sep=""))

print(paste0("dim(cellMETA):",dim(cellMETA)))
print(Sys.time())

##############
# EXP META

library(writexl)

print("**********")
print("writing excel file")
print(Sys.time())

flds=c("title",
"summary",
"dataset_type",
"annotation_source",
"annotation_release_number",
"geo_accession",
"contact_email",
"contact_institute",
"contact_name",
"sample_taxid",
"sample_organism",
"platform_id",
"instrument_model",
"library_selection",
"library_source",
"library_strategy",
"units",
"pubmed_id",
"tags")

vals=c(
nemoMETA.title,#paste(nemoMETA.title,"[",sampLAB,"]",sep=""),
paste(nemoMETA.summary," : [",sampLAB,"]"," This dataset contains ",(dim(exprs)[2])-1," columns/cells/samples, and ",dim(exprs)[1]," rows/genes/features. ",sep=""),
nemoMETA.dataset_type,
nemoMETA.annotation_source,
nemoMETA.annotation_release_number,
nemoMETA.geo_accession,
nemoMETA.contact_email,
nemoMETA.contact_institute,
nemoMETA.contact_name,
nemoMETA.sample_taxid,
nemoMETA.sample_organism,
nemoMETA.platform_id,
nemoMETA.instrument_model,
nemoMETA.library_selection,
nemoMETA.library_source,
nemoMETA.library_strategy,
nemoMETA.units,
nemoMETA.pubmed_id,
nemoMETA.tags
)

meta=data.frame(matrix(nrow=length(flds),ncol=2,data=cbind(flds,vals)))
colnames(meta)=c("field","value")
print(".xlsx file out")
print(paste(baseDIR,sampLAB,"_NeMO_meta.xlsx",sep=""))
write_xlsx(list(metadata=meta),
           path = paste(baseDIR,sampLAB,"_NeMO_meta.xlsx",sep=""),
		   format_headers = F,use_zip64 = T)

# .tar ball
if(TARball.flag){
print("*************************")
print("making tar ball.")
dirNOW=getwd()
setwd(baseDIR)
system(paste0("tar -czvf ",sampLAB,".tar.gz *",sampLAB,"*.tab"))
setwd(dirNOW)
}

} # add on for NeMO

if(!NeMO.flag&TARball.flag){print("You have requested a TARball for NeMO upload, but yo have not requested NeMO datat prcessing, so there is nothing to put in a TARball.")}


print("************************************************************************************************")
print("************************************************************************************************")
print(sampLAB)
print(Sys.time())
print("End.")
print("************************************************************************************************")
print("************************************************************************************************")

}
