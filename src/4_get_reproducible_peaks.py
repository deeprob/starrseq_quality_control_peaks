import os
import argparse
import utils as ut


def call_peaks_from_subsampled_bam(
    input_library_prefix, output_library_prefix,
    bam_dir, peak_call_dir,
    subsampling_fraction, subsampling_rseed, 
    starrpeaker_data_dir, 
    input_library_short, output_library_short
    ):
    # subsample bam file for both input and output
    input_bam_file = os.path.join(bam_dir, input_library_short, f"{input_library_prefix}.bam")
    output_bam_file = os.path.join(bam_dir, output_library_short, f"{output_library_prefix}.bam")
    frac_bam = str(subsampling_fraction).split(".")[1]
    out_dir = os.path.join(peak_call_dir, output_library_short, f"{frac_bam}_{subsampling_rseed}")
    os.makedirs(out_dir, exist_ok=True)
    subsampled_input_bam = os.path.join(out_dir, f"{input_library_prefix}.bam") 
    subsampled_output_bam = os.path.join(out_dir, f"{output_library_prefix}.bam")
    ut.subsample_bam_file(input_bam_file, subsampled_input_bam, frac_bam, subsampling_rseed)
    ut.subsample_bam_file(output_bam_file, subsampled_output_bam, frac_bam, subsampling_rseed)
    # call starrpeaker peak calling functions
    peaks_output_dir = os.path.join(out_dir, "starrpeaker")
    ut.call_starrpeaker_peaks_helper(
        peaks_output_dir,
        subsampled_input_bam, 
        subsampled_output_bam,
        starrpeaker_data_dir,
        )
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("input_lib", type=str, help="The input library name as given in the meta file")
    parser.add_argument("output_lib", type=str, help="The output library name as given in the meta file")
    parser.add_argument("bam_dir", type=str, help="Directory where the bam files are stored")
    parser.add_argument("peak_call_dir", type=str, help="The root dir where called peaks will be stored")
    parser.add_argument("frac", type=float, help="The subsampling fraction of the bam file")
    parser.add_argument("rseed", type=int, help="The random seed of subsampling")
    parser.add_argument("starrpeaker_data_dir", type=str, help="Starrpeaker required files data dir")
    parser.add_argument("-t", "--threads", type=int, help="number of threads to use", default=64) # TODO: include threads in main argument

    cli_args = parser.parse_args()
    in_lib_args = ut.create_args(cli_args.meta_file, cli_args.input_lib)
    out_lib_args = ut.create_args(cli_args.meta_file, cli_args.output_lib)
    call_peaks_from_subsampled_bam(
        in_lib_args.library_prefix,
        out_lib_args.library_prefix,
        cli_args.bam_dir,
        cli_args.peak_call_dir,
        cli_args.frac,
        cli_args.rseed,
        cli_args.starrpeaker_data_dir,   
        in_lib_args.library_short,
        out_lib_args.library_short, 
        )
