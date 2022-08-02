import os
import utils as ut

def main(gtf_file, roi_file, store_file):
    exon_df = ut.parse_gtf_to_get_exons(gtf_file)
    roi_df = ut.read_roi_file(roi_file)
    ut.create_negative_controls(exon_df, roi_df, store_file)
    return


if __name__ == "__main__":
    gtf_file = "/data5/deepro/starrseq/main_library/0_region_selection/data/tss/Refseq_hg38_gene_coordinates.txt"
    roi_file = "/data5/deepro/starrseq/main_library/2_quality_control_lib/data/master/IN/master_filtered.bed"
    store_file = "/data5/deepro/starrseq/main_library/4_quality_control_peaks/data/negative_controls.bed"
    main(gtf_file, roi_file, store_file)
