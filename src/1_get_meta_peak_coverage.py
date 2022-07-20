# make a meta peak coverage file by calculating the read depth of 
# all the enhancer peaks called in each library across all libraries

import argparse
import os
import utils as ut
import itertools
import pandas as pd

def main(
    peak_dir,
    bam_dir,
    store_dir,
    libs_prefix,
    libs_reps,
    libs_short,
    in_lib_short,
    ):

    # get peak files
    peak_files = [ut.get_peak_file(peak_dir, lib_short) for lib_short in libs_short if lib_short!=in_lib_short]
    # store meta peak file
    meta_peak_file = os.path.join(store_dir, "peak_cov", "meta_peak_file.bed")
    ut.store_meta_peak_file(peak_files, meta_peak_file)
    # get bam files
    lib_bam_files = [ut.get_rep_bam_file(bam_dir, lib_short, lib_prefix, lib_reps) for lib_short, lib_prefix, lib_reps in zip(libs_short, libs_prefix, libs_reps)]
    # get coverage dfs
    lib_cov_dfs = [ut.get_replicate_wise_cov_df(meta_peak_file, bam_files) for bam_files in lib_bam_files]    
    # store coverage dfs
    cov_filename = os.path.join(store_dir, "peak_cov", "meta_cov.csv")
    lib_cov_meta_df = pd.concat(lib_cov_dfs, axis=1)
    # repeat lib shorts for the number of replicates present per library
    lib_short_repeat = list(itertools.chain.from_iterable(itertools.repeat(x, len(y.split())) for x,y in zip(libs_short, libs_reps)))
    lib_cov_meta_df.columns = [f"{l}_{c}" for l,c in zip(lib_short_repeat, lib_cov_meta_df.columns)]
    lib_cov_meta_df.to_csv(cov_filename, index=True)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("peak_dir", type=str, help="Directory where the peak files are stored")
    parser.add_argument("bam_dir", type=str, help="Directory where the bam files are stored")
    parser.add_argument("store_dir", type=str, help="Output files dir")
    parser.add_argument("--libs", type=str, help="The library name as given in the meta file", nargs="+", default=["input"])
    parser.add_argument("--in_lib", type=str, help="The library name of the input file", default="input")

    cli_args = parser.parse_args()
    in_args = ut.create_args(cli_args.meta_file,cli_args.in_lib)
    lib_args = [ut.create_args(cli_args.meta_file,lib) for lib in cli_args.libs]

    main(
        cli_args.peak_dir,
        cli_args.bam_dir,
        cli_args.store_dir,
        [lib.library_prefix for lib in lib_args],
        [lib.library_reps for lib in lib_args],
        [lib.library_short for lib in lib_args],
        in_args.library_short
    )
