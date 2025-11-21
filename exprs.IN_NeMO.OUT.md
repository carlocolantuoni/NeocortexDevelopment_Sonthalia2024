# Standardize and reformat the Multi-Omics data using NeMO Function

  This document explains how to prepare inputs and run the **NeMO.TRIM function** (`Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R`). The output of the (`Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R`)
  can then be directly uploaded onto (`NeMOanalytics.org`)

---
## 1. Requirement
Install and load required libraries:

```r
library(devtools)
install_github("CHuanSite/SJD")

library(dplyr)
library(Seurat)
library(SJD)
```

## 2. Prepare the inputs

### 2.1 Expression matrix
The expression matrix stores the processed count data, formatted as genes × cells/samples (i.e., rows represent gene names, columns represent sample or cell identifiers).

Unlike the `Counts.IN.SeuratDimRedNeMO.OUT.R function`, the `Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R` function does not include a built-in normalization argument for raw datasets. Therefore, it is essential to ensure that your dataset has already been fully processed before extracting and providing the expression matrix as input.

Here, we demonstrate two commonly encountered examples for scRNA-seq and BulkRNA-seq, respectively. 
- scRNA-seq: If your dataset has been processed in Seurat and saved as a Seurat object, the normalized count matrix stored in the assay’s "data" slot is already formatted as genes × cells (depending on what version of Seurat you used).
- Bulk RNA-seq: Processed bulk RNA-seq data are typically stored as genes × samples. As long as the dataset has been properly normalized, it can be used directly as the expression matrix input.

```r
# To extract the normalized count matrix from a fully processed Seurat object
exprs_all_cell_types <- as.matrix(seuratobj@assays$RNA@data)
# Check the dimensions of the expression matrix 
dim(exprs_all_cell_types) 
```

### 2.2 cellMETA data
The cellMETA data stores all sample- and cell-level metadata associated with your dataset. This can include experimental group labels, genotypes, additional annotations, cell type labels (for scRNA-seq), or spatial coordinates (for spatial profiling). In short, anything you know about your samples or cells can be captured here.

Each column in cellMETA represents a variable (categorical or continuous) that can later be used to customize visualizations in the NeMO platform. Crucially, cellMETA should also include dimensional reduction coordinates (e.g., PCA, UMAP, t-SNE), so that these embeddings can be carried over and used for corresponding visualization within NeMO.

For scRNAseq data that has been processed with Seurat, the metadata can be extracted from the metadata slot, and the dimensionality reduction coordinates can be easily accessed under the reduction slots. The first two columns of embedding coordinates can then be concatenated together with the original metadata. 

```r
cell_metadata_all_cell_types <- cbind(seuratobj@meta.data,
                                      seuratobj@reductions$pca@cell.embeddings[,1:2],
                                      seuratobj@reductions$umap@cell.embeddings[,1:2],
                                      seuratobj@reductions$tsne@cell.embeddings[,1:2]) 
```

### 2.3 geneMETA data

The geneMETA data is the returned dataset of the `getMatch` function. The getMatch function from the SJD package maps gene identifiers between formats or species. Given a list of gene names, it returns a metadata table with two columns: the corresponding gene symbols and Ensembl IDs. The input species and identifier type (e.g., symbol, ensembl) must be specified to ensure accurate mapping. 

```r
gene_metadata_all_cell_types <- getMatch(rownames(logNormalized_10k_all_cell_types), 
                                         inSpecies = 'mouse', 
                                         inType = 'symbol', 
                                         newSpecies = 'mouse')
```

## 3. Run the function

`Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R` prepares processed expression data and associated metadata for upload to the NeMO Analytics platform. It takes in a gene expression matrix, cell/sample metadata, and gene metadata, and then packages these inputs into NeMO-compatible outputs. The function generates both an .xlsx metadata file and a .tar.gz data archive, which can be uploaded to NeMO Platform. 

The next section will explain each of the arguments in the `Counts.IN.SeuratDimRedNeMO.OUT.TRIM.R` function

- sampLAB: A short tag used to label output files.
- exprs: The processed expression matrix (genes × samples), with normalized values (not raw counts).
- cellMETA: A table of sample/cell metadata
- geneMETA: A table of gene metadata, typically including gene symbols and Ensembl IDs.
- geneMETA.symbol.column: Column name in geneMETA containing gene symbols 
- geneMETA.ensembl.column: Column name in geneMETA containing Ensembl gene IDs 
- baseDIR: Path to the output directory where files will be saved.
- dropDUPna.flag: Logical; if TRUE, duplicate or missing values in metadata are dropped.
- cnvrtDUPna.flag: Logical; if TRUE, duplicate or missing values are converted/handled automatically.
- nemoMETA.title: Title for the dataset (appears in NeMO).
- nemoMETA.summary: Short description or summary of the dataset.
- nemoMETA.dataset_type: Type of dataset (e.g., "scRNA-seq", "BulkRNA-seq", "spatial-transcriptomics").
- nemoMETA.annotation_source: Annotation reference source (e.g., "Ensembl").
- nemoMETA.annotation_release_number: Version number of the annotation release (default = "103").
- nemoMETA.geo_accession: GEO accession number if available.
- nemoMETA.contact_email: Contact email for the dataset.
- nemoMETA.contact_institute: Affiliated institution.
- nemoMETA.contact_name: Contact person’s name.
- nemoMETA.sample_taxid: NCBI Taxonomy ID for the organism (e.g., 9606 = human, 10090 = mouse).
- nemoMETA.sample_organism: Organism name (e.g., "Homo sapiens", "Mus musculus").
- nemoMETA.platform_id: Platform identifier (optional).
- nemoMETA.instrument_model: Sequencing instrument model (optional).
- nemoMETA.library_selection: Library selection method (e.g., "cDNA", "polyA").
- nemoMETA.library_source: Source material for library prep (e.g., "transcriptomic").
- nemoMETA.library_strategy: Strategy used for library prep (e.g., "RNA-Seq").
- nemoMETA.units: Measurement units (e.g., "logCPM", "logTPM").
- nemoMETA.pubmed_id: PubMed ID for an associated publication. 
- nemoMETA.tags: Keywords or tags describing the dataset (comma-separated). This will appear on NeMO as a searchable keyword. 
- NeMO.flag: If TRUE, outputs are formatted for NeMO upload.
- TARball.flag: If TRUE, generates a .tar.gz archive of the expression data.

```r
CntsIN.seurNeMOoutTRIM(
    sampLAB ="Green_Nature_2024_microglia", 
    exprs = exprs_all_cell_type,           
    cellMETA= cellMETA_all_cell_types,      
    geneMETA = geneMETA_all_cell_types,    
    geneMETA.symbol.column="genes",
    geneMETA.ensembl.column="Gene.stable.ID.human",      
    baseDIR=/path/to/file/
    cnvrtDUPna.flag=TRUE,
    nemoMETA.title="Title for the dataset,
    nemoMETA.summary="short description for the study and dataset, preferably including the study design, number of control and        experimental groups, number of cell types, number of cells, etc. ",
    nemoMETA.dataset_type="scRNA-seq",
    nemoMETA.annotation_source="Ensembl",
    nemoMETA.annotation_release_number="103",
    nemoMETA.geo_accession="NA",
    nemoMETA.contact_email="",
    nemoMETA.contact_institute="",
    nemoMETA.contact_name="",
    nemoMETA.sample_taxid=9606,# 10090 for mouse
    nemoMETA.sample_organism="Homo sapiens",# Mus musculus
    nemoMETA.platform_id="",
    nemoMETA.instrument_model="Illumina NovaSeq 6000",
    nemoMETA.library_selection="",
    nemoMETA.library_source="",
    nemoMETA.library_strategy="",
    nemoMETA.units="LogCPM",
    nemoMETA.pubmed_id="39198642",
    nemoMETA.tags="Human, scRNA-seq, AD, cognitive function, ageing")

```





