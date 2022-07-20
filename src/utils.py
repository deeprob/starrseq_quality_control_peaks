import os
import json
from argparse import Namespace
import multiprocessing as mp
import pandas as pd
import pysam
import pybedtools

pybedtools.helpers.set_tempdir("/data5/deepro/tmp")

###############################
# read meta file; create args #
###############################

def create_args(meta_file, lib_name):
    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)
        
    args = Namespace(
        # from metadata file
        library_prefix = meta_dict[lib_name]["prefix"],
        library_reps = meta_dict[lib_name]["replicates"],
        library_pair= meta_dict[lib_name]["read_pairs"],
        library_umi = meta_dict[lib_name]["umi"],
        library_suffix = meta_dict[lib_name]["suffix"],
        library_short = meta_dict[lib_name]["shortform"],
        reference_genome = meta_dict["genome"]["ref_fasta"],
        reference_genome_twobit = meta_dict["genome"]["ref_twobit"],
        roi_file = meta_dict["roi"]["sorted"]
    )

    return args

####################
# filename parsing #
####################

def get_peak_file(peak_dir, lib_short, method="starrpeaker"):
    peak_file_dict = {
        "starrpeaker": "peaks.peak.final.bed",
        "cradle": "CRADLE_peaks",
        "macs2": "NA_peaks.narrowPeak"
    }
    return os.path.join(peak_dir, lib_short, method, peak_file_dict[method])

def get_rep_bam_file(bam_dir, lib_short, lib_prefix, lib_reps):
    bam_files = [os.path.join(bam_dir, lib_short, f"{'_'.join([lib_prefix, rep])}.bam") for rep in lib_reps.split()]
    return bam_files

def get_rep_cov_filename(store_dir, lib_short, filebase):
    return os.path.join(store_dir, "peak_cov", lib_short, f"{filebase}.csv")

#####################
# peak bed coverage #
#####################

def read_bed(bed_file):
    return pybedtools.BedTool(bed_file)

def get_bam_file_reads(bam_file):
    bam_counts = pysam.view("-c", bam_file).strip()
    return int(bam_counts)

def get_replicate_norm_cov(peak_file, bam_file):
    peak_bed = pybedtools.BedTool(peak_file)
    bam_bed = pybedtools.BedTool(bam_file)
    bam_counts = get_bam_file_reads(bam_file)
    cov_bed = peak_bed.coverage(bam_bed)
    cov_df = cov_bed.to_dataframe(disable_auto_names=True, header=None)
    cov_df = cov_df.set_index([0,1,2])
    cov_reads = cov_df.iloc[:, -4]
    cov_reads_rpm = (cov_reads*1e6)/bam_counts
    return cov_reads_rpm

def get_replicate_wise_cov_df(peak_file, bam_files):
    pool_iter = [(peak_file, bf) for bf in bam_files]
    rep_cov_sers = multi_args_pool_job(get_replicate_norm_cov, pool_iter)
    rep_cov_dfs = pd.concat(rep_cov_sers, axis=1)
    rep_cov_dfs.columns = [f"R{i}" for i in range(1, rep_cov_dfs.shape[1] + 1)]
    return rep_cov_dfs

##################
# meta peak file #
##################

def store_meta_peak_file(peak_files, store_path):
    peak_beds = [pybedtools.BedTool(pf) for pf in peak_files]
    meta_peak_bed = pybedtools.BedTool("", from_string=True).cat(*peak_beds, postmerge=False).sort()
    meta_peak_bed.moveto(store_path)
    return

################
# multiprocess #
################

def multi_args_pool_job(func, pool_iter):
    pool = mp.Pool(len(pool_iter))
    res = pool.starmap(func, pool_iter)
    pool.close()
    pool.join()
    return res

