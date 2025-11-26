# h5ad_NeMO: Memory adaptive converter for NeMO format

This module provides a set of Python and R scripts to inspect large `h5ad` files, standardize their metadata, and convert the expression matrix into NeMO compatible text files. The main script can handle `h5ad` files of virtually any size by adapting to the available memory, and it also supports cell filtering and downsampling based on metadata.

See more detailed guide in `h5ad_NeMO_usage.pptx`.

---

## Repository contents

- `install.py`  
  One time installation script that installs all required Python dependencies into the user local environment.

- `h5ad_info.py`  
  Utility script for inspecting the structure of an `h5ad` file, including `X`, `obs`, and `var`.

- `run.py`  
  User facing wrapper script. You only need to edit parameters in this file and then run it from the terminal.

- `h5ad_NeMO_VariableMemory.py`  
  Core implementation of the memory adaptive conversion logic.

- `gene_converter.R`  
  R helper script that converts gene symbols to Ensembl IDs (Version 130).

- `Rfilter_process.py`  
  Additional processing and filtering utilities.

- `h5Seurat_to_h5ad.R`  
  Example script for converting an `h5Seurat` object to `h5ad`.

---

## Setup with `install.py`
```bash
python install.py
```

---

## Inspect .h5ad file with `h5ad_info.py` 
Before running the converter, it is recommended to inspect the input h5ad file and identify: - the data type and format of the **X** matrix - the column in **obs** that stores cell identifiers - the columns in **var** that store Ensembl IDs and gene symbols
```bash
python h5ad_info.py -i path/to/your_file.h5ad
```
In the output: 
- In the **"X MATRIX"** section, confirm whether the data are raw counts or normalized values. 
- In the **"OBSERVATIONS (obs)"** section, find the column that contains cell IDs (for example _index). 
- In the **"VARIABLES (var)"** section, find the column names that correspond to Ensembl IDs and gene symbols (for example ensembl and symbol).
---

## Run the converter with `run.py`
`run.py` is designed so that users do not have to modify the core implementation. You only need to edit parameters in a single place, then execute:
```bash
python run.py
```
Below is a detailed description of all commonly used parameters.

| Parameter                          | Type  | Required | Default                | Description |
|------------------------------------|-------|----------|------------------------|-------------|
| `h5adfile`                         | str   | Yes      | none                   | Path to the input `.h5ad` file. |
| `samp_lab`                         | str   | Yes      | none                   | Output file prefix used in all NeMO text files. |
| `exprs`                            | str   | Yes       | `"X"`                  | Expression matrix to use, typically `"X"` (cells × genes) or a layer name. |
| `cell_meta`                        | str   | Yes       | `"obs"`                | Name of the cell metadata table in the `h5ad` object. |
| `cell_meta_observation_column`     | str   | Yes       | `"_index"`             | Column in `cell_meta` containing cell IDs; used in `_COLmeta.tab` and `_DataMTX.tab`. |
| `cell_filter`                      | str   | No       | `""`                   | R-style logical filter expression for cell selection. |
| `downsample_cell_filter`           | list  | No       | `[]`                   | Column names used for stratified downsampling. |
| `downsample_ratio`                 | float | No       | `1.0`                  | Fraction of cells to retain (0.0–1.0). `1.0` means no downsampling. |
| `gene_meta`                        | str   | Yes       | `"var"`                | Name of the gene metadata table in the `h5ad` object. |
| `gene_meta_symbol_column`          | str   | Yes       | `"Gene"`               | Column in `gene_meta` containing gene symbols. |
| `gene_meta_ensembl_column`         | str   | Yes       | `"Accession"`          | Column in `gene_meta` containing Ensembl gene IDs. |
| `symbol_to_ensembl_rscript`        | str   | No       | `"gene_converter.R"`   | R script for mapping gene symbols to Ensembl IDs. |
| `base_dir`                         | str   | No       | `"./"`                 | Output directory for all NeMO files and intermediates. |
| `drop_dup_na_flag`                 | bool  | No       | `False`                | If `True`, drop genes with duplicated or NA identifiers. |
| `cnvrt_dup_na_flag`                | bool  | No       | `True`                 | If `True`, convert duplicated or NA gene IDs to unique placeholder names. |
| `nemo_meta_title`                  | str   | No       | `""`                   | Title in NeMO metadata. |
| `nemo_meta_summary`                | str   | No       | `""`                   | Summary or abstract of the dataset. |
| `nemo_meta_dataset_type`           | str   | No       | `"scRNA-seq"`          | Dataset type description. |
| `nemo_meta_annotation_source`      | str   | No       | `"Ensembl"`            | Source of gene annotation. |
| `nemo_meta_annotation_release_number` | str | No    | `"103"`                | Annotation release version. |
| `nemo_meta_geo_accession`          | str   | No       | `""`                   | GEO accession ID. |
| `nemo_meta_contact_email`          | str   | No       | `""`                   | Contact email. |
| `nemo_meta_contact_institute`      | str   | No       | `""`                   | Contact institution. |
| `nemo_meta_contact_name`           | str   | No       | `""`                   | Contact name. |
| `nemo_meta_sample_taxid`           | int   | No       | `9606`                 | NCBI taxonomy ID (`9606` human, `10090` mouse). |
| `nemo_meta_sample_organism`        | str   | No       | `"Homo sapiens"`       | Scientific name of the sample organism. |
| `nemo_meta_platform_id`            | str   | No       | `""`                   | Sequencing platform ID. |
| `nemo_meta_instrument_model`       | str   | No       | `""`                   | Sequencing instrument model. |
| `nemo_meta_library_selection`      | str   | No       | `""`                   | Library selection method. |
| `nemo_meta_library_source`         | str   | No       | `""`                   | Library source type. |
| `nemo_meta_library_strategy`       | str   | No       | `""`                   | Library strategy (e.g., RNA-Seq, snRNA-seq). |
| `nemo_meta_units`                  | str   | No       | `""`                   | Units of expression values (e.g., counts, CPM, TPM). |
| `nemo_meta_pubmed_id`              | str   | No       | `""`                   | PubMed ID. |
| `nemo_meta_tags`                   | str   | No       | `""` | Comma-separated dataset tags. |
| `nemo_flag`                        | bool  | No       | `True`                 | Whether to generate NeMO metadata output. |
| `tar_ball_flag`                    | bool  | No       | `True`                 | Whether to compress the final output into a tarball. |
| `memory_size`                      | int   | Yes       | `5`                    | Available memory in GB; determines full-load vs chunked mode. |

**p.s.** Please adjust `memory_size` according to your actual computing environment.
This parameter determines whether the converter uses full-load mode or chunked processing, and setting it too low may lead to slow performance, while setting it too high may exceed available system memory.

