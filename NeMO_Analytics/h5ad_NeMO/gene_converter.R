# Multi-species gene ID conversion
# usage: Rscript gene_converter.R <input_json> <output_json> <input_type> <species>

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]
input_type <- args[3]
species <- args[4]

library(biomaRt)
library(jsonlite)

species_mapping = data.frame(
        species = c("homo sapiens","mus musculus", "caenorhabditis elegans", "drosphila melanagaster",
                    "danio rerio", "xenopus tropicalis", "gallus gallus", "ratus norvegicus", "cavia porcellus",
                    "melanochromis auratus", "oryctolagus cuniculus", "sus scrofa domesticus",
                    "ovis aries", "bos taurus", "canis lupus familiaris",
                    "felis catus", "macaca mulatta", "pan paniscus", "pan troglodytes"),
        ensembl.nms = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl",
                        "celegans_gene_ensembl", "dmelanogaster_gene_ensembl",
                        "drerio_gene_ensembl", "xtropicalis_gene_ensembl", "ggallus_gene_ensembl", "rnorvegicus_gene_ensembl",
                        "cporcellus_gene_ensembl", "mauratus_gene_ensembl",
                        "ocuniculus_gene_ensembl", "sscrofa_gene_ensembl",
                        "oaries_gene_ensembl", "btaurus_gene_ensembl", "clfamiliaris_gene_ensembl",
                        "fcatus_gene_ensembl", "mmulatta_gene_ensembl", "ppaniscus_gene_ensembl",
                        "ptroglodytes_gene_ensembl"),
        genesymbol.attr = c("hgnc_symbol",
                            "mgi_symbol", "external_gene_name", "external_gene_name",
                            "hgnc_symbol", "hgnc_symbol", "hgnc_symbol", "external_gene_name",
                            "hgnc_symbol", "mgi_symbol", "hgnc_symbol", "hgnc_symbol",
                            "hgnc_symbol", "hgnc_symbol", "hgnc_symbol", "hgnc_symbol",
                            "hgnc_symbol", "hgnc_symbol", "hgnc_symbol"))
rownames(species_mapping) <- species_mapping$species

species <- tolower(species)
if (!species %in% rownames(species_mapping)) {
    stop(paste("Error: Species", species, "not found in mapping table"))
}
dataset_name = species_mapping[species, ]$ensembl.nms
mart <- useEnsembl(biomart = "ensembl", dataset = dataset_name, 
                   host = "https://may2021.archive.ensembl.org")
cat(paste("Start executing", input_type, "Convert...\n"))

if (input_type == "symbol") {
    filter = species_mapping[species, ]$genesymbol.attr
    output = "ensembl_gene_id"
}else if (input_type == "ensembl") {
    filter = "ensembl_gene_id"
    output = species_mapping[species, ]$genesymbol.attr
} else {
    stop("Error: input_type must be 'symbol' or 'ensembl'")
}

genes <- fromJSON(input_path)
unique_genes <- unique(genes)
result <- getBM(attributes = c(filter, output), 
             filters = filter, values = unique_genes, mart = mart)

output_genes <- rep("n/a", length(genes))

for (i in 1:length(genes)) {
    gene <- genes[i]
    match_rows <- which(result[, filter] == gene)
    
    if (length(match_rows) > 0) {
        matched_results <- result[match_rows, , drop = FALSE]
        
        if (nrow(matched_results) == 1) {
            # Single match - directly assign
            conversion_result <- matched_results[1, output]
            if (!is.na(conversion_result) && conversion_result != "") {
                output_genes[i] <- as.character(conversion_result)
            }
        } else {
            # Multiple matches - need to select one
            valid_results <- matched_results[!is.na(matched_results[, output]) & 
                                           matched_results[, output] != "", , drop = FALSE]
            
            if (nrow(valid_results) > 0) {
                # Sort valid results and take the first one (alphabetically smallest)
                output_values <- valid_results[, output]
                sorted_values <- sort(output_values, na.last = TRUE)
                selected_value <- sorted_values[1]
                
                if (!is.na(selected_value) && selected_value != "") {
                    output_genes[i] <- as.character(selected_value)
                    cat(paste("Info: Multiple matches found for", gene, 
                             "- selected:", selected_value, "\n"))
                } else {
                    cat(paste("Warning: Multiple matches found for", gene, 
                             "but no valid results\n"))
                }
            } else {
                cat(paste("Warning: Multiple matches found for", gene, 
                         "but all results are empty/NA\n"))
            }
        }
    }
    # If no matches found, output_genes[i] remains "n/a"
}

failed_count <- sum(output_genes == "n/a")
success_count <- length(genes) - failed_count
cat(paste("Successfully converted:", success_count, "genes\n"))
cat(paste("Failed to convert:", failed_count, "genes (marked as n/a)\n"))

write_json(output_genes, output_path, pretty = TRUE, auto_unbox = TRUE)
cat("Conversion completed")