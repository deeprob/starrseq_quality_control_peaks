# store negative control file coverage df
# dir structure: data/peak_cov/lib_short/{input|output}.csv

import argparse
import os
import utils as ut


def main(
    negative_control_file,
    bam_dir,
    store_dir,
    in_prefix,
    in_reps,
    in_short,
    lib_prefix,
    lib_reps,
    lib_short,
    ):

    # get input bam files
    in_bam_files = ut.get_rep_bam_file(bam_dir, in_short, in_prefix, in_reps)
    # get output bam files
    lib_bam_files = ut.get_rep_bam_file(bam_dir, lib_short, lib_prefix, lib_reps)
    # get coverage dfs
    in_cov_df = ut.get_replicate_wise_cov_df(negative_control_file, in_bam_files)
    lib_cov_df = ut.get_replicate_wise_cov_df(negative_control_file, lib_bam_files)      
    # store coverage dfs
    in_cov_filename, out_cov_filename = ut.get_neg_control_cov_filename(store_dir, lib_short, "input"), ut.get_neg_control_cov_filename(store_dir, lib_short, "output")
    os.makedirs(os.path.dirname(in_cov_filename), exist_ok=True)
    in_cov_df.to_csv(in_cov_filename, index=True)
    lib_cov_df.to_csv(out_cov_filename, index=True)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib", type=str, help="The library name as given in the meta file")
    parser.add_argument("negative_control_file", type=str, help="The filepath to the negative control file")
    parser.add_argument("bam_dir", type=str, help="Directory where the bam files are stored")
    parser.add_argument("store_dir", type=str, help="Output files dir")
    parser.add_argument("--in_lib", type=str, help="The input library name as given in the meta file", default="input")

    cli_args = parser.parse_args()
    in_args = ut.create_args(cli_args.meta_file, cli_args.in_lib)
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)

    main(
        cli_args.negative_control_file,
        cli_args.bam_dir,
        cli_args.store_dir,
        in_args.library_prefix,
        in_args.library_reps,
        in_args.library_short,
        lib_args.library_prefix, 
        lib_args.library_reps,
        lib_args.library_short
    )
