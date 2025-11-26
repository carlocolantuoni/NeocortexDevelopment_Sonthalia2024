library(Seurat)
library(SeuratDisk)
library(reticulate)

##  Input Path (all h5Seurat files under the input path will be converted) ##
input_dict <- ""

## Output Path (default: the current dict) ##
output_dir <- getwd()


if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 
files <- list.files(path = input_dict, pattern = "\\.h5Seurat$", full.names = TRUE)

for (f in files) {
  cat("Processing h5Seurat file:", f, "\n")

  obj <- LoadH5Seurat(f, assays = "RNA",
                      reductions = FALSE,
                      graphs = FALSE,
                      neighbors = FALSE,
                      images = FALSE,
                      tools = FALSE,
                      misc = FALSE)

  DefaultAssay(obj) <- "RNA"
  obj[["RNA"]]@scale.data <- matrix(numeric(), nrow = 0, ncol = 0)
  obj@reductions <- list()
  obj@graphs <- list()
  obj@neighbors <- list()
  obj@tools <- list()
  obj@commands <- list()
  obj@misc <- list()

  out_name <- sub("\\.h5Seurat$", ".h5ad", basename(f))
  temp_h5seurat <- file.path(output_dir, paste0("temp_", basename(f)))
  SaveH5Seurat(obj, filename = temp_h5seurat, overwrite = TRUE)

  Convert(temp_h5seurat, dest = "h5ad", overwrite = TRUE)

  temp_h5ad <- sub("\\.h5Seurat$", ".h5ad", temp_h5seurat)
  out_path <- file.path(output_dir, out_name)
  file.rename(temp_h5ad, out_path)

  if (file.exists(temp_h5seurat)) {
    unlink(temp_h5seurat)
  }
}

