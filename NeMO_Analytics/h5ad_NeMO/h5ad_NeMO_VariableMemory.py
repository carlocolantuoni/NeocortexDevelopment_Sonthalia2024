import os
import re
import numpy as np
import subprocess
import time
import gc
from tqdm import tqdm
from openpyxl import Workbook
import json
import h5py
import scanpy as sc
from scipy.sparse import csr_matrix
import inspect
from sklearn.model_selection import train_test_split


def counts_in_h5ad_nemo_out_trim(
        h5adfile,
        samp_lab: str,  # output file prefix
        exprs: str = 'X',  # cells * genes
        cell_meta: str = 'obs',
        cell_meta_observation_column: str = "_index",  # Column name in cell_meta for observation IDs
                                                    # used as the first column in _COLmeta.tab and the colname in _DataMTX.tab
        cell_filter: str="",     # R-style filter expression for cell selection
        downsample_cell_filter: list =[],   # List of column names for stratified downsampling
        downsample_ratio: float=1.0,    # Fraction of cells to retain during downsampling (0.0-1.0). 
                                        # 1.0 means no downsampling
        
        gene_meta: str = 'var',
        gene_meta_symbol_column: str = "Gene",  # Column name in gene_meta for gene symbols
        gene_meta_ensembl_column: str = "Accession",  # Column name in gene_meta for ensembl IDs
        symbol_to_ensembl_rscript: str = "gene_converter.R", # TSV file, containing a mapping table of Ensembl gene IDs and gene symbols
        base_dir: str = "/dcs05/carlo/legacy-dcl01-ccolantu/data/Explr/",  # output dir
        drop_dup_na_flag: bool = False,
        cnvrt_dup_na_flag: bool = True,
        nemo_meta_title: str = "",
        nemo_meta_summary: str = "",
        nemo_meta_dataset_type: str = "scRNA-seq",
        nemo_meta_annotation_source: str = "Ensembl",
        nemo_meta_annotation_release_number: str = "103",
        nemo_meta_geo_accession: str = "",
        nemo_meta_contact_email: str = "",
        nemo_meta_contact_institute: str = "",
        nemo_meta_contact_name: str = "",
        nemo_meta_sample_taxid: int = 9606,  # 10090 for mouse
        nemo_meta_sample_organism: str = "Homo sapiens",  # Mus musculus
        nemo_meta_platform_id: str = "",
        nemo_meta_instrument_model: str = "",
        nemo_meta_library_selection: str = "",
        nemo_meta_library_source: str = "",
        nemo_meta_library_strategy: str = "",
        nemo_meta_units: str = "",
        nemo_meta_pubmed_id: str = "",
        nemo_meta_tags: str = "tissue, cells, experiment",
        nemo_flag: bool = True,
        tar_ball_flag: bool = True,
        memory_size: int = 5
) -> None:
    """
    Convert expression matrix (potentially from h5ad) to NeMO format files.
    """

    ### Parameter Validation ###
    frame = inspect.currentframe()
    local_vars = frame.f_locals
    for k, v in local_vars.items():
        if isinstance(v, str):
            local_vars[k] = v.strip()

    # No spaces in the samp_lab to avoid issues with tar command
    samp_lab = samp_lab.replace(" ", "")

    # Ensure base_dir ends with a slash
    os.makedirs(base_dir, exist_ok=True)
    if not base_dir.endswith('/'):
        base_dir += '/'

    # Create output directory if it doesn't exist
    os.makedirs(base_dir, exist_ok=True)

    print("*" * 80)
    print("*" * 80)
    print("Begin.")
    print(samp_lab)
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    print("*" * 80)
    print("*" * 80)

    # Check if title will be too long for NeMO
    if nemo_flag and len(f"{nemo_meta_title}[{samp_lab}]") > 255:
        print("ERROR: Length of final pasted NeMO title is too long (>255 characters). Ending processing.")
        print(f"{nemo_meta_title}[{samp_lab}]")
        return None

    # Validate input parameters
    if drop_dup_na_flag and cnvrt_dup_na_flag:
        print("ERROR: Cannot request both drop_dup_na_flag AND cnvrt_dup_na_flag")
        return None

    if not nemo_flag and tar_ball_flag:
        print(
            "You have requested a TARball for NeMO upload, but you have not requested NeMO data processing, so there is nothing to put in a TARball.")
        return None

    if not nemo_flag:
        print("NeMO output not requested. Exiting.")
        return None

    # Name temporary files
    exprs_file = f"{base_dir}{samp_lab}_DataMTX.tab"
    gene_meta_file = f"{base_dir}{samp_lab}_ROWmeta.tab"
    cell_meta_file = f"{base_dir}{samp_lab}_COLmeta.tab"
    temp_exprs_file = f"{base_dir}{samp_lab}_DataMTX_initial.tab"
    temp_cell_meta_file = f"{base_dir}{samp_lab}_COLmeta_initial.tab"
    filtered_h5ad_file = f"{base_dir}{samp_lab}_filtered.h5ad"

    #####################
    # Filtering
    if len(cell_filter)!=0:
        def h5ad_filter(h5ad_file,r_expr,output_path):
            from Rfilter_process import translate_r_filter
            query_expr = translate_r_filter(r_expr)
            adata = sc.read_h5ad(h5ad_file)
            adata.raw = None
            filtered_obs = adata.obs.query(query_expr)
            adata_filtered = adata[filtered_obs.index].copy()
            del adata
            del filtered_obs
            gc.collect()
            print("Start saving the filtered h5ad file")
            adata_filtered.write_h5ad(output_path)
            print("Complete Save")
            return output_path
        h5adfile = h5ad_filter(h5adfile,cell_filter.strip(),filtered_h5ad_file)
    
    # Downsampling
    if downsample_ratio < 1.0:
        def h5ad_downsample(h5ad_file,downsample_cell_filter,downsample_ratio,output_path):
            adata = sc.read_h5ad(h5ad_file)
            adata.raw = None
            if len(downsample_cell_filter)>0:
                missing_cols = [col for col in downsample_cell_filter if col not in adata.obs.columns]
                if missing_cols:
                    raise ValueError(f"The following columns don't exist: {missing_cols}")
                if len(downsample_cell_filter) == 1:
                    stratify_col = downsample_cell_filter[0]
                    sampled_obs = adata.obs.groupby(stratify_col, group_keys=False).apply(
                        lambda group: group.sample(
                            n=max(1, int(len(group) * downsample_ratio)),
                            random_state=42
                        ) if len(group) > 0 else group
                    )
                else:
                    stratify_labels = adata.obs[downsample_cell_filter].apply(
                        lambda row: '|'.join(row.astype(str)), axis=1
                    )
                    indices = np.arange(len(adata.obs))
                    _, sampled_indices = train_test_split(
                        indices,
                        test_size=downsample_ratio,  # test_size是我们要保留的比例
                        stratify=stratify_labels,
                        random_state=42
                    )
                    sampled_obs = adata.obs.iloc[sampled_indices]
            else:
                sampled_obs = adata.obs.sample(frac=downsample_ratio,random_state=42)
            adata_downsampled = adata[sampled_obs.index].copy()
            adata_downsampled.write_h5ad(output_path)
            return output_path
        h5adfile = h5ad_downsample(h5adfile,downsample_cell_filter,
                                  downsample_ratio,filtered_h5ad_file)

    #####################
    # Add on for NeMO
    print("**********")
    print("NeMO output requested.")
    print(samp_lab)
    print(time.strftime("%Y-%m-%d %H:%M:%S"))

    def get_h5ad_path(string_list):
        path_list = []
        for string in string_list:
            path_parts = string.split('/')
            path_parts = [f'"{item.strip()}"' for item in path_parts]
            path_list.append('f['+']['.join(path_parts)+']')
        return path_list

    def identify_underlying_structure(group, group_string):
        """Identify the underlying structure of an HDF5 data set"""
        attri_dict = {'group': group, 'row': None, 'col': None, 'is_sparse': False, 'format': None}
        if isinstance(group, h5py.Group):
            if all(k in group for k in ['data', 'indices', 'indptr']):
                attri_dict['is_sparse'] = True
                indptr_shape = group['indptr'].shape
                # Determine if it's CSR or CSC
                matrix_shape = group.attrs.get('shape')
                if matrix_shape is not None:
                    n_rows, n_cols = matrix_shape
                    n_ptr = indptr_shape[0]
                    attri_dict['row'] = n_rows
                    attri_dict['col'] = n_cols
                    if n_ptr == n_rows + 1:
                        attri_dict['format'] = 'CSR'
                    elif n_ptr == n_cols + 1:
                        attri_dict['format'] = 'CSC'
                    else:
                        attri_dict['format'] = 'CSR'
                        print(f"\033[93mWarning\033[0m: {group_string} is a sparse matrix,\
                         but the encoding format is unknown. Try to process it according to the default CSR")
            else:
                subgroups = list(group.keys())
                num_datasets = 0
                single_dimensional = True
                for k in subgroups:
                    if isinstance(group[k], h5py.Dataset):
                        num_datasets += 1
                        if len(group[k].shape) != 1:
                            single_dimensional = False
                if len(subgroups) == 2 and num_datasets == 2 and single_dimensional:
                    attri_dict['col'] = 0
                    if all(k in subgroups for k in ['categories', 'codes']):
                        attri_dict['format'] = subgroups
                        attri_dict['row'] = len(group['codes'])
                    else:
                        categorical = False
                        if all((0 <= item < len(group[subgroups[0]])) for item in group[subgroups[1]]):
                            categorical = True
                            attri_dict['format'] = ['categories', 'codes']
                            attri_dict['row'] = len(group[subgroups[0]])
                        else:
                            if all((0 <= item < len(group[subgroups[1]])) for item in group[subgroups[0]]):
                                categorical = True
                                attri_dict['format'] = ['codes', 'categories']
                                attri_dict['row'] = len(group[subgroups[1]])
                        if not categorical:
                            print(f"\033[93mWarning\033[0m: {group_string} contains two datasets but\
                                                         the categories format is not recognized, so skip this group")
                            return None
                else:
                    print(f"\033[93mWarning\033[0m: {group_string} contains multiple datasets, so skip this group.")
                    return None

        elif isinstance(group, h5py.Dataset):
            if len(group.shape) == 1:
                attri_dict['row'] = len(group)
                attri_dict['col'] = 0
            elif len(group.shape) == 2:
                attri_dict['row'], attri_dict['col'] = group.shape
            else:
                print(f"\033[93mWarning\033[0m: {group_string} unkown structure, so skip this group.")
                return None
        return attri_dict

    def decode_if_bytes(dataset):
        """Decode byte data into string format"""
        data = dataset[:]
        if isinstance(data, np.ndarray):
            data = data.tolist()
        while isinstance(data[0], bytes):
            data = [x.decode('utf-8').strip("b'\"") for x in data]
        if not isinstance(data[0], str):
            data = [str(x) for x in data]
        data = np.array(data)
        return data

    def Dataset_Processing(group, format):
        """Handling datasets of different formats (categorical or normal)"""
        if isinstance(group, h5py.Group):
            subgroups = list(group.keys())
            cat_idx, code_idx = (0, 1) if format == ['categories', 'codes'] else (1, 0)
            categories = decode_if_bytes(group[subgroups[cat_idx]])
            codes = group[subgroups[code_idx]][:]
            data = categories[codes]
        else:
            data = decode_if_bytes(group)
        return data

    def append_row_to_tab(data, temp_file):
        """Append a row of data to the tab file"""
        if os.path.exists(temp_file):
            with open(temp_file, 'a') as f_in:
                f_in.write('\t'.join(data) + '\n')
        else:
            with open(temp_file, 'w') as f_in:
                f_in.write('\t'.join(data) + '\n')

    def transpose_tab(input_file, output_file, memory_size):
        """Transpose tab file"""
        num_row = 1
        with open(input_file, 'r') as f_in:
            num_col = len(f_in.readline().strip().split('\t'))
            for _ in f_in:
                num_row += 1
        file_size = os.path.getsize(input_file) /(1024**3)
        buffer_size_mb = min(512, max(1, int(memory_size * 0.05 * 1024)))
        if 0.7*memory_size > file_size:
            data_array = np.genfromtxt(input_file, delimiter='\t', dtype=str)
            transposed = data_array.T
            with open(output_file, 'a', buffering=buffer_size_mb * 1024 * 1024) as f_out:
                for row in transposed:
                    f_out.write('\t'.join(row) + '\n')
            del data_array, transposed
            gc.collect()
        else:
            chunk = max(1, int(0.5*memory_size*num_row//file_size))
            total_chunks = (num_col - 1) // chunk + 1
            with tqdm(total=num_col, desc="Transpose the matrix", unit="Col") as pbar:
                for col_start in range(0, num_col, chunk):
                    col_end = min(col_start + chunk, num_col)
                    current_chunk_size = col_end - col_start
                    chunk_array = np.empty((num_row, current_chunk_size), dtype=str)
                    with open(input_file, 'r') as f_in:
                        for i, line in enumerate(f_in):
                            parts = line.strip().split('\t')
                            chunk_array[i, :] = parts[col_start:col_end]
                    chunk_transposed = chunk_array.T
                    with open(output_file, 'a', buffering=buffer_size_mb * 1024 * 1024) as f_out:
                        for data in chunk_transposed:
                            f_out.write('\t'.join(data) + '\n')
                    pbar.update(current_chunk_size)
                    pbar.set_postfix({"Chnk": f"{col_start // chunk + 1}/{total_chunks}"})
            del chunk_array, chunk_transposed
            gc.collect()

    def MergeFiles_ByCol(temp_file, merged_indices, output_file, buffer_size):
        """Merge multiple temporary files by column"""
        file_handles = [open(temp_file, 'r')]
        file_handles += [open(f"{temp_file}_merged_{index}", 'r') for index in merged_indices]
        with open(temp_file, 'r') as f:
            total_lines = sum(1 for _ in f)
        with open(output_file, 'w') as out_file:
            line_count = 0
            with tqdm(total=total_lines, desc="Merging the Files") as pbar:
                while True:
                    batch_lines = []
                    for _ in range(buffer_size):
                        lines = []
                        eof = False
                        # Read a line from each file
                        for fh in file_handles:
                            line = fh.readline()
                            if not line:
                                eof = True
                                break
                            lines.append(line.strip())
                        if eof:
                            break
                        batch_lines.append('\t'.join(lines))
                    # Batch Write
                    if batch_lines:
                        out_file.write('\n'.join(batch_lines) + '\n')
                        line_count += len(batch_lines)
                        pbar.update(len(batch_lines))
                    if eof or len(batch_lines) < buffer_size:
                        break
                    # Perform garbage collection every 10 batches
                    if (line_count // buffer_size) % 10 == 0:
                        gc.collect()
        # Close all the files
        for fh in file_handles:
            fh.close()
        os.remove(temp_file)
        for index in merged_indices:
            os.remove(f"{temp_file}_merged_{index}")

    def get_dtype_size(group):
        """Get the data type size of each component of the sparse matrix for memory estimation"""
        dataset_list = ['data','indices','indptr']
        length = [len(group[item]) for item in dataset_list]
        try:
            type = [str(group[item].dtype) for item in dataset_list]
            pattern = re.compile(r'\d+')
            size = [int(pattern.search(item).group())//8 for item in type]
            return sum(np.array(length)*np.array(size))/(1024*1024*1024), size[0]
        except:
            if len(group['data'])>2.1 * 10**9:
                size = 8
            else:
                size = 4
            return sum(np.array(length)*np.array([4,size,size]))/(1024*1024*1024),4

    def CSR_to_tab(group, final_file, temp_file, ensembl_ids, memory_size):
        """Convert CSR format sparse matrix to tab file"""
        cell_num, gene_num = group.attrs.get('shape')
        CSR_memory,value_byte = get_dtype_size(group)
        buffer_size_mb = min(512, max(1, int(memory_size * 0.05 * 1024)))
        str_vectorized = np.vectorize(str)
        if memory_size*0.35 > CSR_memory:
            data = group['data'][:]
            indices = group['indices'][:]
            indptr = group['indptr'][:]
            csr = csr_matrix((data, indices, indptr), shape=(cell_num, gene_num))
            print("Convert CSR to CSC format...")
            csc = csr.tocsc()
            del csr, data, indices, indptr
            gc.collect()

            with open(final_file, 'a', buffering = buffer_size_mb * 1024 * 1024) as f_out:
                with tqdm(total=gene_num, desc="Writing gene expression data", unit="gene") as pbar:
                    for gene in range(gene_num):
                        start_idx = csc.indptr[gene]
                        end_idx = csc.indptr[gene + 1]
                        gene_data = np.full(cell_num, '0', dtype=str)
                        if start_idx < end_idx:
                            row_indices = csc.indices[start_idx:end_idx]
                            values = csc.data[start_idx:end_idx]
                            gene_data[row_indices] = str_vectorized(values)
                        f_out.write(ensembl_ids[gene] + '\t' + '\t'.join(gene_data) + '\n')
                        pbar.update(1)

        else:
            with open(temp_file, 'w') as f_exprs:
                for _ in range(gene_num):
                    f_exprs.write(ensembl_ids[_] + '\n')
            temp_indices = []
            merged_indices = []
            value_memory = 2*cell_num*gene_num*value_byte/(1024*1024*1024)
            cell_chunk = max(1, int(0.1*memory_size*cell_num/(CSR_memory+value_memory)))
            total_cell_chunks = (cell_num - 1) // cell_chunk + 1
            gene_chunk = max(1, int(0.1*memory_size*gene_num/(CSR_memory+value_memory)))
            merge_every = min(50, 1024*1024//(len(temp_file) + 55))
            with tqdm(total=cell_num, desc="Processing CSR Matrix", unit="Cell") as pbar:
                for index, cell_start in enumerate(range(0, cell_num, cell_chunk)):
                    cell_end = min(cell_start + cell_chunk, cell_num)
                    current_chunk_size = cell_end - cell_start
                    # Extract CSR data in chunks
                    indptr = group['indptr'][cell_start:cell_end + 1]
                    start_index = indptr[0]
                    end_index = indptr[-1]
                    data = group['data'][start_index:end_index]
                    indices = group['indices'][start_index:end_index]
                    indptr = indptr - start_index
                    # Transpose and save
                    sub_csr = csr_matrix((data, indices, indptr), shape=(current_chunk_size, gene_num))
                    chunk_file = f"{temp_file}_{index}"
                    temp_indices.append(index)
                    np.savetxt(chunk_file, sub_csr.toarray().T, delimiter='\t', fmt='%s')
                    # Update progress bar
                    pbar.update(current_chunk_size)
                    pbar.set_postfix({"Chunk": f"{index + 1}/{total_cell_chunks}"})
                    if len(temp_indices) >= merge_every or cell_end == cell_num:
                        merge_index = len(merged_indices)
                        merged_file = f"{temp_file}_merged_{merge_index}"
                        merged_indices.append(merge_index)
                        temp_files_paths = [f"{temp_file}_{i}" for i in temp_indices]
                        try:
                            if os.name == 'posix':
                                files_str = ' '.join(temp_files_paths)
                                cmd = f"paste {files_str} > {merged_file}"
                                subprocess.run(cmd, shell=True, check=True)
                            elif os.name == 'nt':
                                for temp_f in temp_files_paths:
                                    with open(temp_f, 'rb') as infile, open(temp_file, 'ab', buffering=buffer_size_mb  * 1024 * 1024) as outfile:
                                        while True:
                                            chunk = infile.read(buffer_size_mb  * 1024 * 1024)
                                            if not chunk:
                                                break
                                            outfile.write(chunk)
                            else:
                                raise ('Unsupported system')
                        except Exception as e:
                            print(f"Error during file merge: {e}")
                            raise
                        # Clean up temporary files
                        for temp_f in temp_files_paths:
                            try:
                                os.remove(temp_f)
                            except Exception as e:
                                print(f"Warning: Could not remove {temp_f}: {e}")
                        temp_indices = []
            del sub_csr, data, indices, indptr
            gc.collect()
            intermediate_file = temp_file + '.final'
            MergeFiles_ByCol(temp_file, merged_indices, intermediate_file, gene_chunk)
            with open(intermediate_file, 'rb') as source, open(final_file, 'ab', buffering=buffer_size_mb * 1024 * 1024) as dest:
                while True:
                    chunk = source.read(buffer_size_mb* 1024 * 1024)
                    if not chunk:
                        break
                    dest.write(chunk)
            os.remove(intermediate_file)

    def CSC_to_tab(group, output_file, ensembl_ids, memory_size):
        """Convert CSC format sparse matrix to tab file"""
        cell_num, gene_num = group.attrs.get('shape')
        buffer_size_mb = min(512, max(1, int(memory_size * 0.05 * 1024)))
        str_vectorized = np.vectorize(str)
        with open(output_file, 'a', buffering=buffer_size_mb * 1024 * 1024) as f_out:
            with tqdm(total=gene_num, desc="Writing gene expression data", unit="gene") as pbar:
                for gene in range(gene_num):
                    start_idx = group['indptr'][gene]
                    end_idx = group['indptr'][gene + 1]
                    gene_data = np.full(cell_num, '0', dtype=str)
                    if start_idx < end_idx:
                        row_indices = group['indices'][start_idx:end_idx]
                        values = group['data'][start_idx:end_idx]
                        gene_data[row_indices] = str_vectorized(values)
                    f_out.write(ensembl_ids[gene] + '\t' + '\t'.join(gene_data) + '\n')
                    pbar.update(1)
        del gene_data, row_indices, values
        gc.collect()

    def DenseMatrix_to_tab(dataset, output_file, ensembl_ids, memory_size):
        """Convert dense matrix to tab file"""
        cell_num, gene_num = dataset.shape
        value_byte = dataset.dtype.itemsize
        value_memory = cell_num*gene_num*(value_byte*4+100)/(1024*1024*1024)
        buffer_size_mb = min(512, max(1, int(memory_size * 0.05 * 1024)))
        if 0.7*memory_size > value_memory:
            value = dataset[:].T
            with open(output_file, 'a', buffering=buffer_size_mb * 1024 * 1024) as f_out:
                with tqdm(total=gene_num, desc="Writing gene expression data", unit="gene") as pbar:
                    for i, row in enumerate(value):
                        f_out.write(ensembl_ids[i] + '\t' + '\t'.join(str(x) for x in row) + '\n')
                        pbar.update(1)
            del value
            gc.collect()
        else:
            batch_size = max(1, int(0.5*memory_size*cell_num/value_memory))
            with open(output_file, 'a', buffering=buffer_size_mb * 1024 * 1024) as f_out:
                with tqdm(total=gene_num, desc="Writing gene expression data", unit="gene") as pbar:
                    for batch_start in range(0, gene_num, batch_size):
                        batch_end = min(batch_start + batch_size, gene_num)
                        batch_data = dataset[:, batch_start:batch_end].T
                        for i, row in enumerate(batch_data):
                            gene_idx = batch_start + i
                            f_out.write(ensembl_ids[gene_idx] + '\t' +
                                        '\t'.join(str(x) for x in row) + '\n')
                            pbar.update(1)
                    del batch_data
                    gc.collect()

    def dup_na(ensembl, symbol):
        """"Identify gene IDs and symbols for duplications and deletions"""
        seen_ensembl = set()
        ensembl_dup_na = []
        na_values = ['n/a', 'None', 'NA', 'N/A', 'NaN', 'nan', 'null', '']
        for i, ens in enumerate(ensembl):
            if ens in na_values or ens.isspace():
                ensembl_dup_na.append(i)
            elif ens in seen_ensembl:
                ensembl_dup_na.append(i)
            else:
                seen_ensembl.add(ens)

        seen_symbols = set()
        symbol_dup_na = []
        for i, sym in enumerate(symbol):
            if sym in na_values or sym.isspace():
                symbol_dup_na.append(i)
            elif sym in seen_symbols:
                symbol_dup_na.append(i)
            else:
                seen_symbols.add(sym)

        return ensembl_dup_na, symbol_dup_na

    def delete_rows_from_file(filename, rows_to_delete):
        """Delete the specified line from the file (from the end to the front)"""
        batch_size = min(1000, (1024*1024-len(filename))/(2*len(rows_to_delete)))
        adjusted_rows = [r + 2 for r in rows_to_delete]
        adjusted_rows.sort(reverse=True)
        for i in range(0, len(adjusted_rows), batch_size):
            batch = adjusted_rows[i:i + batch_size]
            sed_expr = ';'.join(f"{row}d" for row in batch)
            cmd = f"sed -i '{sed_expr}' {filename}"
            subprocess.run(cmd, shell=True, check=True)

    # ==================== Main data processing flow ==================== #
    ##### Raw data loading ####
    with h5py.File(h5adfile, 'r') as f:
        input_list = [exprs, cell_meta, cell_meta + '/' + cell_meta_observation_column,
                       gene_meta + '/' + gene_meta_symbol_column, gene_meta + '/' + gene_meta_ensembl_column]
        record_list = get_h5ad_path(input_list)
        record = [None for k in record_list]
        ### Create the reference ###
        subgroups = list(f[cell_meta].keys())
        record[1] = {'group': [], 'row': [], 'col': [], 'is_sparse': [], 'format': [], 'subgroup': []}
        for k in subgroups:
            cellmeta_record = f'{record_list[1]}["{k}"]'
            subrecord = identify_underlying_structure(eval(cellmeta_record), cellmeta_record)
            if subrecord is not None:
                if k == cell_meta_observation_column:
                    record[2] = subrecord
                else:
                    record[1]['subgroup'].append(k)
                    for key, val in zip(record[1].keys(), subrecord.values()):
                        record[1][key].append(val)

        record[0] = identify_underlying_structure(eval(record_list[0]), record_list[0])
        record[3] = identify_underlying_structure(eval(record_list[3]), record_list[3])

        if len(gene_meta_ensembl_column)!=0:
            cnvrt_symbol = False
            record[4] = identify_underlying_structure(eval(record_list[4]), record_list[4])
            ## Basic Test ##
            for i in range(len(record)):
                if record[i] is None:
                    print(f"\033[91mError\033[0m: Fail in identifying attributes of {record_list[i]}")
                    return None
            # Test for gene_meta
            gene_nums = record[3]['row']
            if record[3]['row'] != record[4]['row']:
                print(f"gene symbols: {record[3]['row']}")
                print(f"ensembl IDs: {record[4]['row']}")
                print(
                    f"\033[91mError\033[0m: The number of rows in the gene symbols column and the ensembl IDs column are not equal")
                return None
        else:
            cnvrt_symbol = True
            ## Basic Test ##
            for i in range(len(record)-1):
                if record[i] is None:
                    print(f"\033[91mError\033[0m: Fail in identifying attributes of {record_list[i]}")
                    return None
            # Test for gene_meta
            gene_nums = record[3]['row']

        # Test for cell_meta
        cell_nums_list = record[1]['row']
        cell_nums = cell_nums_list[0]
        for nums in cell_nums_list:
            if nums is not None and nums != cell_nums:
                print(f"\033[91mError\033[0m: The number of rows in the cell meta data is not equal")
                return None
        if cell_nums != record[2]['row']:
            print(f"cell num: {cell_nums}")
            print(f"observation: {record[2]['row']}")
            print(f"\033[91mError\033[0m: The number of rows in the cell meta data is not equal")
            return None
        # Test for exprs
        pass_test = True
        if cell_nums != record[0]['row']:
            pass_test = False
            print(f"cells: {cell_nums}")
            print(f"expr_row: {record[0]['row']}")
            print(f"\033[91mError\033[0m: The number of rows and cells in the expression matrix are not equal")
        if record[3]['row'] != record[0]['col']:
            pass_test = False
            print(f"genes: {gene_nums}")
            print(f"expr_col: {record[0]['col']}")
            print(
                f"\033[91mError\033[0m: The number of columns in the expression matrix is not equal to the number of genes")
        if not pass_test:
            return None

        ### Generate col&row names ###
        observation = Dataset_Processing(record[2]['group'], record[2]['format'])
        symbol = Dataset_Processing(record[3]['group'], record[3]['format'])
        # convert symbol to ensembl
        if cnvrt_symbol:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            rscript_path = os.path.join(script_dir, symbol_to_ensembl_rscript)
            if os.path.exists(rscript_path):
                symbol_path = f"{base_dir}{samp_lab}_symbol.json"
                ensembl_path = f"{base_dir}{samp_lab}_ensembl.json"
                with open(symbol_path, 'w') as f_symbol:
                    json.dump(symbol.tolist(), f_symbol)
                command=f"""
                module load R && \
                srun --mem=10G Rscript {rscript_path} "{symbol_path}" "{ensembl_path}" "symbol" "{nemo_meta_sample_organism}"
                """
                convert_result = subprocess.run(command, shell=True, executable="/bin/bash",
                                        capture_output=True, text=True)
                if convert_result.stdout:
                    print(convert_result.stdout)
                if convert_result.stderr:
                    print(convert_result.stderr)
                if os.path.exists(symbol_path):
                    os.remove(symbol_path)
                if convert_result.returncode == 0:
                    with open(ensembl_path,'r') as f_ensembl:
                        ensembl = np.array(json.load(f_ensembl))
                    os.remove(ensembl_path)
                else:
                    print(f"\033[91mError\033[0m: The symbol_to_ensembl convertion failed")
                    return None
            else:
                print(f"\033[91mError\033[0m: The symbol_to_ensembl Rscript doesn't exist")
                return None
        else:
            ensembl = Dataset_Processing(record[4]['group'], record[4]['format'])

        if cnvrt_dup_na_flag or drop_dup_na_flag:
            ensembl_dup_na, symbol_dup_na = dup_na(ensembl, symbol)
            if cnvrt_dup_na_flag:
                if len(ensembl_dup_na) > 0:
                    print(f"Convert {len(ensembl_dup_na)} Ensembl ID duplicates/NA to dummy IDs")
                    replacement_ensembl = np.array([f"NOensemblIDmapped.{i + 1}" for i in range(len(ensembl_dup_na))])
                    ensembl[ensembl_dup_na] = replacement_ensembl
                if len(symbol_dup_na) > 0:
                    print(f"Convert {len(symbol_dup_na)} Gene symbol duplicates/NA to dummy IDs")
                    replacement_symbol = np.array([f"NOsymbolMapped.{i + 1}" for i in range(len(symbol_dup_na))])
                    symbol[symbol_dup_na] = replacement_symbol

        ### Generate temporary files ###
        print("**********")
        print("geneMETA file")
        with open(exprs_file,'w') as f_exprs, open(gene_meta_file, 'w') as f_gene, open(cell_meta_file, 'w') as f_cell:
            # gene_meta
            f_gene.write('gene\tgene_symbol\n')
            for _ in range(gene_nums):
                f_gene.write(ensembl[_]+'\t'+symbol[_]+'\n')
            # headers of cell_meta and exprs
            f_cell.write('observations\t'+'\t'.join(map(str,record[1]['subgroup']))+'\n')
            f_exprs.write('\t'+'\t'.join(map(str,observation))+'\n')
        print(time.strftime("%Y-%m-%d %H:%M:%S"))

        # cell_meta
        print("**********")
        print("cellMETA file")
        append_row_to_tab(observation, temp_cell_meta_file)
        for _ in range(len(record[1]['subgroup'])):
            cell_subdata = Dataset_Processing(record[1]['group'][_],record[1]['format'][_])
            append_row_to_tab(cell_subdata, temp_cell_meta_file)
        transpose_tab(temp_cell_meta_file,cell_meta_file,memory_size)
        os.remove(temp_cell_meta_file)
        del cell_subdata
        print(time.strftime("%Y-%m-%d %H:%M:%S"))

        # exprs
        print("**********")
        print("exprs file")
        if record[0]['is_sparse']:
            if record[0]['format'] == 'CSR':
                CSR_to_tab(record[0]['group'],exprs_file, temp_exprs_file, ensembl, memory_size)
            else:
                CSC_to_tab(record[0]['group'], exprs_file, ensembl, memory_size)
        else:
            DenseMatrix_to_tab(record[0]['group'], exprs_file, ensembl, memory_size)
        print(time.strftime("%Y-%m-%d %H:%M:%S"))

        if drop_dup_na_flag:
            all_dup_na = set(ensembl_dup_na + symbol_dup_na)
            delete_rows_from_file(gene_meta_file, all_dup_na)
            delete_rows_from_file(exprs_file, all_dup_na)
            print(f"Delete {len(all_dup_na)} duplicates/NA rows")

    # Create metadata for excel file
    fields = [
        "title",
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
        "tags"
    ]

    vals = [
        nemo_meta_title,
        f"{nemo_meta_summary} : [{samp_lab}] This dataset contains {cell_nums} columns/cells/samples, and {gene_nums} rows/genes/features.",
        nemo_meta_dataset_type,
        nemo_meta_annotation_source,
        nemo_meta_annotation_release_number,
        nemo_meta_geo_accession,
        nemo_meta_contact_email,
        nemo_meta_contact_institute,
        nemo_meta_contact_name,
        nemo_meta_sample_taxid,
        nemo_meta_sample_organism,
        nemo_meta_platform_id,
        nemo_meta_instrument_model,
        nemo_meta_library_selection,
        nemo_meta_library_source,
        nemo_meta_library_strategy,
        nemo_meta_units,
        nemo_meta_pubmed_id,
        nemo_meta_tags
    ]

    # Create Excel workbook using openpyxl (instead of writexl in R)
    wb = Workbook()
    ws = wb.active
    ws.title = "metadata"

    # Add header
    ws.append(["field", "value"])

    # Add data
    for i, (field, val) in enumerate(zip(fields, vals)):
        ws.append([field, val])

    # Save Excel file
    print(".xlsx file")
    xlsx_file = f"{base_dir}{samp_lab}_NeMO_meta.xlsx"
    wb.save(f"{xlsx_file}")

    # Create tar.gz file
    if tar_ball_flag:
        print("*************************")
        print("making tar ball.")
        files_to_pack = [exprs_file, gene_meta_file, cell_meta_file]
        existing_files = [f for f in files_to_pack if os.path.exists(f)]

        if existing_files:
            try:
                output_file = f"{base_dir}{samp_lab}.tar.gz"
                cmd = ["tar", "-czf", output_file] + existing_files
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    print(f"Created {samp_lab}.tar.gz with {len(existing_files)} files")
                else:
                    print(f"Error creating tar file: {result.stderr}")
            except Exception as e:
                print(f"Error: {e}")
        else:
            print("No files found to archive")

    print("*" * 80)
    print("*" * 80)
    print(samp_lab)
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    print("End.")
    print("*" * 80)
    print("*" * 80)


if __name__ == "__main__":
    # use example
    counts_in_h5ad_nemo_out_trim(
        h5adfile='/dcs05/carlo/legacy-dcl01-ccolantu/data/Explr/Green_Nature_2024_h5ad/inhibitory.h5ad',
        samp_lab="", # output file prefix
        exprs="X",  # "X" or "layers/counts"
        cell_meta='obs',
        cell_meta_observation_column="",  # Column name in cell_meta for observation IDs
        cell_filter="",     # Such as "cell_type %in% c('B cell', 'T cell', 'NK cell')"
        downsample_cell_filter=[],   # Such as ["sample","cell_type"]
        downsample_ratio=0.1,    # (0,1)
        gene_meta='var',
        gene_meta_symbol_column="",  # Column name in gene_meta for gene symbols
        gene_meta_ensembl_column="",  # Column name in gene_meta for ensembl IDs. If "" (empty string), automatically match according to gene symbol
        base_dir="/dcs05/carlo/legacy-dcl01-ccolantu/data/Explr/",  # output dir
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
