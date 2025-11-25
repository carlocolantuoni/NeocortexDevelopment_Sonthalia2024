# NeocortexDevelopment_Sonthalia2024

Repository for analyses performed in "NeMO Analytics: A Compendium of Transcriptomic Data for the Exploration of Neocortical Development" (https://doi.org/10.1101/2024.02.26.581612).

10min video summary of our manuscript: https://youtu.be/IE3itucbig4

Interactive, public data exploration resources can be found at: https://nemoanalytics.org/landing/neocortex/

Data aquisition from NeMO Anlaytics and Joint decomposition performed for analysis in Figures 1, 2 & 3: SJDcode_Figure1n2n3.R

Data aquisition from NeMO Anlaytics and Joint decomposition performed for analysis in Figures 4, 5, & 6: SJDcode_Figure4n5n6.R

Gene loadings from the 3 decompositions in the paper: SonthaliaEtAl_SupplTable1_p7CtxDev_GeneLoadings.txt; SonthaliaEtAl_SupplTable2_p40CtxDev_GeneLoadings.txt; SonthaliaEtAl_SupplTable3_p20CtxLayers_GeneLoadings.txt

R code used for performing gene set enrichment analysis: GSEAcode_SonthaliaEtAl.R

Gene lists used in the GSEA were obtained from: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/ and CtxDiseaseGeneListsPLUS.RData from PMID: 40634286

GeneCart2idCols.R is a function used in the SJD code (SJDcode_Figure1n2n3.R & SJDcode_Figure4n5n6.R) to create files for the import of joint decomposition results into NeMO Analytics as gene lists for projection.

Conda environment and associated packages used in the pesudotime analysis can be recreated using the Seurat5_R4_3_0.yml file (https://github.com/carlocolantuoni/NeocortexDevelopment_Sonthalia2024/blob/main/Seurat5_R4_3_0.yml) and the command "conda env create -f Seurat5_R4_3_0.yml".


Installation instructions and code for joint matrix decomposition using the SJD package (https://www.biorxiv.org/content/10.1101/2022.11.07.515489v1) can be found here: https://github.com/CHuanSite/SJD

Installation instructions and code for projection analysis using projectR (https://academic.oup.com/bioinformatics/article/36/11/3592/5804979) can be found here: https://www.bioconductor.org/packages/release/bioc/html/projectR.html

R function for getting data into a format that can be uploaded to NeMO Analytics: exprs.IN_NeMO.OUT.R

python function for getting .h5ad data into a format that can be uploaded to NeMO Analytics: h5ad_NeMO

youtube playlist for NeMO Analytics tutorials: https://www.youtube.com/playlist?list=PLbeFhDgIkXxe0G1hzEB_vLSP-5qqesNQj


Additional help can be found on the github site for gEAR, the platform on which NeMO Analytics is built (https://github.com/IGS/gEAR/wiki) and on NeMO Analytics webisite itself:
https://nemoanalytics.org/upload_dataset.html; 
https://nemoanalytics.org/dataset_explorer.html; & 
https://nemoanalytics.org/gene_list_manager.html

For more general info on dataset collections and computational tools go to: https://www.carlocolantuoni.org/
