import sys
sys.path.append('/dcs05/carlo/legacy-dcl01-ccolantu/data/Rfunctions/h5ad_NeMO')
from h5ad_NeMO_VariableMemory import counts_in_h5ad_nemo_out_trim

counts_in_h5ad_nemo_out_trim(
    h5adfile='/dcs05/carlo/legacy-dcl01-ccolantu/data/Explr/Green_Nature_2024_cellular_trajectory_ageing_AD/scRNAseq/complete_metadata_h5ad', # h5ad_path
    base_dir="./",  # output dir
    samp_lab="Green_Nature_2024_cellular_trajectory_ageing_AD", # output file prefix
    exprs="X",  # "X" or "layers/counts"
    cell_meta='obs',
    cell_meta_observation_column="",  # Column name in cell_meta for observation IDs
    cell_filter="",     # empty or "cell_type %in% c('B cell', 'T cell', 'NK cell')"
    downsample_cell_filter=[],   # empty or ["sample","cell_type"]
    downsample_ratio=1.0,    # 1.0 or (0.0,1.0)
    gene_meta='var',
    gene_meta_symbol_column="",  # Column name in gene_meta for gene symbols
    gene_meta_ensembl_column="",  # Column name in gene_meta for ensembl IDs. If "" (empty string), automatically match according to gene symbol
    drop_dup_na_flag=False,
    cnvrt_dup_na_flag=True,
    nemo_meta_title="",
    nemo_meta_summary="""
    
    """,
    nemo_meta_dataset_type="scRNA-seq",
    nemo_meta_annotation_source="Ensembl",
    nemo_meta_annotation_release_number="103",
    nemo_meta_geo_accession="",
    nemo_meta_contact_email="",
    nemo_meta_contact_institute="",
    nemo_meta_contact_name="",
    nemo_meta_sample_taxid=9606,  # 10090 for mouse
    nemo_meta_sample_organism="Homo sapiens",  # Mus musculus
    nemo_meta_platform_id="",
    nemo_meta_instrument_model="",
    nemo_meta_library_selection="",
    nemo_meta_library_source="",
    nemo_meta_library_strategy="",
    nemo_meta_units="=",
    nemo_meta_pubmed_id="=",
    nemo_meta_tags="=",
    nemo_flag=True,
    tar_ball_flag=True,
    memory_size = 100
)