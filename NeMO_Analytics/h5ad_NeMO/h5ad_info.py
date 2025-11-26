#!/usr/bin/env python3

import h5py
import os
import numpy as np
import argparse

def decode_if_bytes(data):
    if isinstance(data, np.ndarray):
        data = data.tolist()
    while isinstance(data, list) and data and isinstance(data[0], bytes):
        data = [x.decode('utf-8').strip("b'\"") for x in data]
    return data

def explore_h5ad(filename):
    """Explore the basic structure and attributes of an h5ad file"""
    print(f"{'=' * 60}")
    print(f"{'ANALYZING FILE:':^60}")
    print(f"{filename:^60}")
    print(f"{'=' * 60}")

    file_size = os.path.getsize(filename) / (1024 * 1024)  # MB
    print(f"File size: {file_size:.2f} MB")

    with h5py.File(filename, 'r') as f:
        # Check top-level structure
        print(f"\n{'=' * 60}")
        print(f"{'TOP-LEVEL STRUCTURE':^60}")
        print(f"{'=' * 60}")
        top_level_keys = list(f.keys())
        print(f"Top-level keys: {top_level_keys}")

        # Check X matrix structure
        print(f"\n{'=' * 60}")
        print(f"{'X MATRIX':^60}")
        print(f"{'=' * 60}")
        if 'X' in f:
            if isinstance(f['X'], h5py.Group):
                print("X is a group (sparse matrix)")
                x_keys = list(f['X'].keys())
                print(f"Keys in X: {x_keys}")

                # Check matrix format
                if all(k in f['X'] for k in ['data', 'indices', 'indptr']):
                    data_shape = f['X/data'].shape
                    indices_shape = f['X/indices'].shape
                    indptr_shape = f['X/indptr'].shape
                    print(f"data shape: {data_shape}")
                    print(f"indices shape: {indices_shape}")
                    print(f"indptr shape: {indptr_shape}")

                    # Determine if it's CSR or CSC
                    matrix_shape = f['X'].attrs.get('shape')
                    if matrix_shape is not None:
                        print(f"Matrix shape: {matrix_shape}")
                        n_rows, n_cols = matrix_shape
                        n_ptr = indptr_shape[0]

                        if n_ptr == n_rows + 1:
                            print("Format: CSR (Compressed Sparse Row)")
                        elif n_ptr == n_cols + 1:
                            print("Format: CSC (Compressed Sparse Column)")
                        else:
                            print(
                                f"Unknown format: indptr length ({n_ptr}) doesn't match rows+1 ({n_rows + 1}) or cols+1 ({n_cols + 1})")

                    # View data samples - first 10 elements
                    print("\nData samples:")
                    print("First 10 elements of data:", f['X/data'][0:10])
                    print("First 10 elements of indices:", f['X/indices'][0:10])
                    print("First 10 elements of indptr:", f['X/indptr'][0:10])

            elif isinstance(f['X'], h5py.Dataset):
                print("X is a dataset (dense matrix)")
                print(f"X shape: {f['X'].shape}")
                print(f"X type: {f['X'].dtype}")
                print("X data sample (upper-left corner):")

                try:
                    sample = f['X'][0:min(5, f['X'].shape[0]), 0:min(5, f['X'].shape[1])]
                    print(sample)

                    # Print first 10 elements (flattened)
                    flat_sample = f['X'].flatten()[0:10]
                    print("First 10 elements (flattened):", flat_sample)
                except Exception as e:
                    print(f"Cannot sample X data: {e}")
        else:
            print("File doesn't contain X matrix")

        # Check obs and var
        for section in ['obs', 'var']:
            print(f"\n{'=' * 60}")
            if section == 'obs':
                print(f"{'OBSERVATIONS (obs)':^60}")
            else:
                print(f"{'VARIABLES (var)':^60}")
            print(f"{'=' * 60}")

            if section in f:
                section_keys = list(f[section].keys())
                print(f"\n{section} keys: {section_keys}")

                # Get section size information
                num_datasets = 0
                num_groups = 0
                for k in section_keys:
                    if isinstance(f[section][k], h5py.Dataset):
                        num_datasets += 1
                    elif isinstance(f[section][k], h5py.Group):
                        num_groups += 1

                print(f"{section} contains {num_datasets} datasets and {num_groups} groups")

                # Check index
                if '_index' in f[section]:
                    index_type = f[section]['_index'].dtype
                    print(f"{section} index type: {index_type}")
                    try:
                        index_sample = f[section]['_index'][0:10]  # First 10 elements
                        print(f"{section} index sample (first 10): {index_sample}")
                    except Exception as e:
                        print(f"Cannot sample {section} index: {e}")

                # Get number of items (if possible)
                for k in section_keys:
                    if isinstance(f[section][k], h5py.Dataset):
                        try:
                            print(f"\nNumber of {section} items: {len(f[section][k])}")
                            break  # Only need length from one Dataset
                        except:
                            pass

                # Check each key
                for k in section_keys:
                    print(f"\nExamining {section}/{k}")
                    if isinstance(f[section][k], h5py.Dataset):
                        print(f"{section}/{k} is a dataset, shape: {f[section][k].shape}, type: {f[section][k].dtype}")
                        try:
                            sample_data = f[section][k][0:10]  # First 10 elements
                            sample_data = decode_if_bytes(sample_data)
                            print(f"{section}/{k} sample (first 10): {sample_data}")
                        except Exception as e:
                            print(f"Cannot sample {section}/{k} data: {e}")

                    elif isinstance(f[section][k], h5py.Group):
                        group_keys = list(f[section][k].keys())
                        print(f"{section}/{k} is a group with {len(group_keys)} keys: {group_keys}")

                        # Try to get more info about the group
                        for subkey in group_keys[:3]:  # Show only first 3 subkeys
                            if isinstance(f[section][k][subkey], h5py.Dataset):
                                print(
                                    f"  - {subkey} is a dataset, shape: {f[section][k][subkey].shape}, type: {f[section][k][subkey].dtype}")
                                try:
                                    sample_data = f[section][k][subkey][0:10]  # First 10 elements
                                    sample_data = decode_if_bytes(sample_data)
                                    print(f"    Sample (first 10): {sample_data}")
                                except:
                                    pass
                            elif isinstance(f[section][k][subkey], h5py.Group):
                                subgroup_keys = list(f[section][k][subkey].keys())
                                print(f"  - {subkey} is a nested group with keys: {subgroup_keys[:5]}...")
            else:
                print(f"File doesn't contain {section}")

        # Check other common sections
        print(f"\n{'=' * 60}")
        print(f"{'OTHER COMPONENTS':^60}")
        print(f"{'=' * 60}")
        for component in ['layers', 'uns', 'obsm', 'varm', 'obsp', 'varp']:
            if component in f:
                component_keys = list(f[component].keys())
                print(f"\n{component} keys: {component_keys}")

                for k in component_keys[:5]:  # Show only first 5 keys
                    print(f"\nExamining {component}/{k}")
                    if isinstance(f[component][k], h5py.Dataset):
                        print(
                            f"{component}/{k} is a dataset, shape: {f[component][k].shape}, type: {f[component][k].dtype}")
                        try:
                            if len(f[component][k].shape) <= 1:
                                sample_data = f[component][k][0:10]  # First 10 elements
                                sample_data = decode_if_bytes(sample_data)
                                print(f"{component}/{k} sample (first 10): {sample_data}")
                            else:
                                # For multi-dimensional data, try to get a small slice
                                if len(f[component][k].shape) == 2:
                                    sample_data = f[component][k][0:min(3, f[component][k].shape[0]),
                                                  0:min(3, f[component][k].shape[1])]
                                    sample_data = decode_if_bytes(sample_data)
                                    print(f"{component}/{k} sample (3x3 corner):\n{sample_data}")
                                else:
                                    print(
                                        f"{component}/{k} is multi-dimensional data, showing first 10 flattened elements")
                                    flat_sample = f[component][k].flatten()[0:10]
                                    flat_sample = decode_if_bytes(flat_sample)
                                    print(f"First 10 flattened elements: {flat_sample}")
                        except Exception as e:
                            print(f"Cannot sample {component}/{k} data: {e}")

                    elif isinstance(f[component][k], h5py.Group):
                        subkeys = list(f[component][k].keys())
                        print(f"{component}/{k} is a group with {len(subkeys)} keys: {subkeys[:5]}...")

                        # Try to get subkey information
                        for subkey in subkeys[:3]:  # Show only first 3 subkeys
                            if isinstance(f[component][k][subkey], h5py.Dataset):
                                print(f"  - {subkey} is a dataset, shape: {f[component][k][subkey].shape}")
                                try:
                                    sample_data = f[component][k][subkey][0:10]  # First 10 elements
                                    sample_data = decode_if_bytes(sample_data)
                                    print(f"    Sample (first 10): {sample_data}")
                                except Exception as e:
                                    print(f"    Cannot sample data: {e}")
                            elif isinstance(f[component][k][subkey], h5py.Group):
                                subsub_keys = list(f[component][k][subkey].keys())
                                print(f"  - {subkey} is a nested group with {len(subsub_keys)} keys")

                if len(component_keys) > 5:
                    print(f"... and {len(component_keys) - 5} more keys not shown")
            else:
                print(f"File doesn't contain {component}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Explore h5ad file structure and contents')
    parser.add_argument('-i', '--input', required=True, help='Input h5ad file path')
    args = parser.parse_args()
    explore_h5ad(args.input)
