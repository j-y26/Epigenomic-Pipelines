# Given an input H5 file, extract only the gex (gene expression) features
# and save them to a new H5 file.

# Usage: python extract_gex_from_h5.py -i/--input <input_h5_file>
#                                      -o/--output <output_h5_file>

import argparse
import h5py
import numpy as np
from scipy.sparse import csc_matrix

def extract_gex_from_multiome_h5(input_file, output_file):
    """
    Extract only gene expression data from Cell Ranger ARC multiome h5 file
    while preserving the exact Cell Ranger v3 format structure.
    
    Args:
        input_file (str): Path to input multiome h5 file
        output_file (str): Path to output gene expression only h5 file
    """
    
    with h5py.File(input_file, "r") as f_in:
        # Read matrix components
        barcodes = f_in["matrix/barcodes"][:]
        data = f_in["matrix/data"][:]
        indices = f_in["matrix/indices"][:]
        indptr = f_in["matrix/indptr"][:]
        shape = f_in["matrix/shape"][:]

        # Read features
        feature_type = f_in["matrix/features/feature_type"][:]
        feature_ids = f_in["matrix/features/id"][:]
        feature_names = f_in["matrix/features/name"][:]
        feature_genomes = f_in["matrix/features/genome"][:]
        feature_intervals = f_in["matrix/features/interval"][:]

        # Identify gene expression features
        gex_mask = np.array([ft.decode("utf-8") == "Gene Expression" for ft in feature_type])
        gex_indices = np.where(gex_mask)[0]

        # Filter data for gene expression features
        gex_feature_type = feature_type[gex_mask]
        gex_feature_ids = feature_ids[gex_mask]
        gex_feature_names = feature_names[gex_mask]
        gex_feature_genomes = feature_genomes[gex_mask]
        gex_feature_intervals = feature_intervals[gex_mask]

        # Create a sparse CSC matrix
        matrix = csc_matrix((data, indices, indptr), shape=shape)

        # Filter the matrix to keep only gene expression features
        # This will create a new matrix with only the GEX features
        gex_matrix = matrix[gex_indices, :]

        # Extract the new components from the filtered matrix
        new_data = gex_matrix.data.astype(np.int32)
        new_indices = gex_matrix.indices.astype(np.int64)
        new_indptr = gex_matrix.indptr.astype(np.int64)
        new_shape = (gex_indices.size, shape[1])

        # Create a new HDF5 file
        with h5py.File(output_file, "w") as f_out:
            # Create matrix group
            matrix_group = f_out.create_group("matrix")
            matrix_group.create_dataset("barcodes", data=barcodes, compression="gzip", compression_opts=4)
            matrix_group.create_dataset("data", data=new_data, dtype=np.int32, compression="gzip", compression_opts=4)
            matrix_group.create_dataset("indices", data=new_indices, dtype=np.int64, compression="gzip", compression_opts=4)
            matrix_group.create_dataset("indptr", data=new_indptr, dtype=np.int64, compression="gzip", compression_opts=4)
            matrix_group.create_dataset("shape", data=new_shape, dtype=np.int32, compression="gzip", compression_opts=4)

            # Create features group
            features_group = matrix_group.create_group("features")
            features_group.create_dataset("feature_type", data=gex_feature_type, compression="gzip", compression_opts=4)
            features_group.create_dataset("id", data=gex_feature_ids, compression="gzip", compression_opts=4)
            features_group.create_dataset("name", data=gex_feature_names, compression="gzip", compression_opts=4)
            features_group.create_dataset("genome", data=gex_feature_genomes, compression="gzip", compression_opts=4)
            features_group.create_dataset("interval", data=gex_feature_intervals, compression="gzip", compression_opts=4)

            # Copy additional feature metadata if present
            if "_all_tag_keys" in f_in["matrix/features"]:
                features_group.create_dataset("_all_tag_keys", data=f_in["matrix/features/_all_tag_keys"][:], compression="gzip", compression_opts=4)

            # Copy target_sets if present (only for GEX-related target sets, if applicable)
            if "target_sets" in f_in["matrix/features"]:
                target_sets_group = features_group.create_group("target_sets")
                for target_set in f_in["matrix/features/target_sets"]:
                    # Copy target sets only if they correspond to GEX features
                    # This assumes target_sets are associated with GEX; adjust if needed
                    target_sets_group.create_dataset(target_set, data=f_in["matrix/features/target_sets"][target_set][:], compression="gzip", compression_opts=4)

    print(f"GEX-only HDF5 file created at: {output_file}")


def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description="Given an input H5 file, extract only the gex (gene expression) features and save them to a new H5 file.")

    parser.add_argument("-i", "--input", dest="input_h5", help="Input H5 file", default="raw_feature_bc_matrix.h5")
    parser.add_argument("-o", "--output", dest="output_h5", help="Output H5 file", default="raw_feature_bc_matrix_gex.h5")
    
    args = parser.parse_args()

    extract_gex_from_multiome_h5(args.input_h5, args.output_h5)

if __name__ == "__main__":
    main()