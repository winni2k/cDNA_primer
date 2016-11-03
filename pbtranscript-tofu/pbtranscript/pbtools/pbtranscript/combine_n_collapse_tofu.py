__author__ = 'etseng@pacb.com'

import os, sys, glob
from cPickle import *
import pbtools.pbtranscript.tofu_wrap as tofu_wrap

def find_split_dirs(dirname):
    split_dirs = []
    for d in glob.glob(os.path.join(dirname, '*kb_part*')):
        d2 = os.path.join(dirname, d)
        split_dirs.append(d2)
    return split_dirs


def main(args):
    args.root_dir = os.path.abspath(args.root_dir)
    split_dirs = find_split_dirs(args.root_dir)
    print >> sys.stderr, "The following subdirectories are to be combined:"
    for d in split_dirs: print >> sys.stderr, d

    combined_dir = os.path.join(args.root_dir, 'combined')
    if not os.path.exists(combined_dir):
        os.makedirs(combined_dir)
    # (4) combine quivered HQ/LQ results
    hq_filename, lq_filename, hq_pre_dict, lq_pre_dict = \
            tofu_wrap.combine_quiver_results(split_dirs, combined_dir, args.quiver_hq_filename, args.quiver_lq_filename,\
            args.tofu_prefix)
    with open(os.path.join(args.root_dir, 'combined', 'combined.hq_lq_pre_dict.pickle'), 'w') as f:
        dump({'HQ': hq_pre_dict, 'LQ': lq_pre_dict}, f)

    print >> sys.stderr, "Size bins combined and written to", os.path.join(args.root_dir, 'combined')

    if args.skip_gmap:
        print >> sys.stderr, "Skipping GMAP step. DONE."
        return
    # (5) collapse quivered HQ results
    print >> sys.stderr, "Running GMAP alignment then collapse..."
    collapse_prefix_hq = tofu_wrap.run_collapse_sam(hq_filename, args.gmap_db, args.gmap_name, cpus=args.gmap_nproc, max_fuzzy_junction=args.max_fuzzy_junction, dun_merge_5_shorter=True)
    # (6) make abundance
    tofu_wrap.get_abundance(collapse_prefix_hq, hq_pre_dict, collapse_prefix_hq)
    # (7) run filtering & removing subsets in no5merge
    if args.targeted_isoseq:
        tofu_wrap.run_filtering_by_count(collapse_prefix_hq, collapse_prefix_hq+'.min_fl_5', min_count=5)
        tofu_wrap.run_filtering_away_subsets(collapse_prefix_hq+'.min_fl_5', collapse_prefix_hq+'.min_fl_5.filtered', args.max_fuzzy_junction)
    else:
        tofu_wrap.run_filtering_by_count(collapse_prefix_hq, collapse_prefix_hq+'.min_fl_2', min_count=2)
        tofu_wrap.run_filtering_away_subsets(collapse_prefix_hq+'.min_fl_2', collapse_prefix_hq+'.min_fl_2.filtered', args.max_fuzzy_junction)


if __name__ == "__main__":
    import argparse
    import binascii
    parser = argparse.ArgumentParser(prog='tofu_wrap')

    parser.add_argument("--root_dir", required=True, help="Root dir containing the size bins (ex: clusterOut/)")
    parser.add_argument("--quiver_hq_filename", default="all_quivered_hq.100_30_0.99.fastq", help="Quiver HQ filename (default: all_quivered_hq.100_30_0.99.fastq)")
    parser.add_argument("--quiver_lq_filename", default="all_quivered_lq.fastq", help="Quiver LQ filename (default: all_quivered_lq.fastq)")

    parser.add_argument("--targeted_isoseq", default=False, action="store_true", help="tofu was run with --targeted_isoseq (default: off)")
    parser.add_argument("--skip_gmap", default=False, action="store_true", help="Skip GMAP part")
    parser.add_argument("--gmap_nproc", default=12, type=int, help="GMAP CPU (default: 12)")
    parser.add_argument("--gmap_name", default="hg19", help="GMAP DB name (default: hg19)")
    parser.add_argument("--gmap_db", default="/home/UNIXHOME/etseng/share/gmap_db_new/", help="GMAP DB location (default: /home/UNIXHOME/etseng/share/gmap_db_new/)")
    parser.add_argument("--output_seqid_prefix", type=str, default=None, help="Output seqid prefix. If not given, a random ID is generated")
    parser.add_argument("--max_fuzzy_junction", default=5, type=int, help="Max fuzzy junction (default: 5 bp)")

    args = parser.parse_args()

    args.tofu_prefix = binascii.b2a_hex(os.urandom(3)) if args.output_seqid_prefix is None else args.output_seqid_prefix


    main(args)
