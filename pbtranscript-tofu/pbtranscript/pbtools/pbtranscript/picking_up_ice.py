#!/usr/bin/env python
import os, sys
from cPickle import *
from pbcore.io.FastaIO import FastaReader
import pbtools.pbtranscript.ice.ProbModel as pm
import pbtools.pbtranscript.ice.IceIterative as ice
from pbtools.pbtranscript.ice.IceUtils import ice_fa2fq

def ensure_pickle_goodness(pickle_filename, root_dir, fasta_files_to_add=None):
    """
    Old versions of IceIterative.write_pickle is missing some key/values.
    Add if needed.
    Return a good pickle object
    """
    a = load(open(pickle_filename))
    if a['fasta_filename'] != os.path.abspath(os.path.join(root_dir,'current.fasta')):
        raise Exception, "The pickle file {0} indicates that current.fasta is not being used. ICE likely did not finish to a point that could be picked up.".format(pickle_filename)

    a['newids'] = check_n_fix_newids(a)

    if fasta_files_to_add is not None:
        for file in fasta_files_to_add.split(','):
            if not os.path.exists(file):
                raise Exception, "{0} is not a valid fasta file to add!".format(file)
            if file in a['fasta_filenames_to_add']:
                print >> sys.stderr, "{0} is already in to-add list. Ignore.".format(file)
            a['fasta_filenames_to_add'].append(file)
    if 'root_dir' not in a:
        print >> sys.stderr, "Pickle {0} missing some key-values. Fixing it.".format(pickle_filename)
        a['root_dir'] = root_dir
        a['all_fasta_filename'] = a['all_fasta_fiilename']
        a['qv_prob_threshold'] = 0.03
        with open(pickle_filename + '.fixed', 'w') as f:
            dump(a, f)
        print >> sys.stderr, "Fixed pickle written to {0}.fixed".format(pickle_filename)
        return a, pickle_filename + '.fixed'
    else:
        # newid might have been fixed, STILL output pickle writing anyway
        with open(pickle_filename, 'w') as f:
            dump(a, f)
        return a, pickle_filename

def check_n_fix_newids(icec_obj):
    newids = icec_obj['newids']

    if len(newids) == 0:
        print >> sys.stderr, "newids is empty (probably a finalized run). set it."
        for k, v in icec_obj['d'].iteritems():
            if len(v) != 1:
                newids.add(k)
        print >> sys.stderr, "added {0} seqs to newids".format(len(newids))
    return newids

def make_current_fasta(icec_obj, flnc_filename, root_dir):
    """
    current fasta will consists of all ids

    however --- if this was a already finished run and we are adding more input,
        then newids is empty, in this case we set newids = everything that
        has no affiliation or more than one affiliated cluster in d
    """
    with open(os.path.join(root_dir, 'current.fasta'), 'w') as f:
        for r in FastaReader(flnc_filename):
                f.write(">{0}\n{1}\n".format(r.name, r.sequence))

def pickup_icec_job(pickle_filename, ccs_fofn, flnc_filename, fasta_files_to_add, root_dir):
    icec_obj, icec_pickle_filename = ensure_pickle_goodness(pickle_filename, root_dir, fasta_files_to_add)
    make_current_fasta(icec_obj, flnc_filename, root_dir)
    print >> sys.stderr, "Reading QV information...."
    # first need to convert to fastq
    ice_fa2fq('current.fasta', ccs_fofn, 'current.fastq')
    probqv = pm.ProbFromFastq(os.path.join(root_dir,'current.fastq'))
    icec = ice.IceIterative.from_pickle(icec_pickle_filename, probqv)
    # first must RE-RUN gcon to get all the proper refs
    icec.changes = set()
    icec.refs = {}
    icec.ccs_fofn = ccs_fofn
    icec.all_fasta_filename = flnc_filename
    todo = icec.uc.keys()
    print >> sys.stderr, "Re-run gcon for proper refs...."
    icec.run_gcon_parallel(todo)
    print >> sys.stderr, "Re-calculating cluster prob, just to be safe...."
    icec.calc_cluster_prob(True)
    print >> sys.stderr, "Sanity checking now...."
    icec.sanity_check_uc_refs()
    icec.ensure_probQV_newid_consistency()
    print >> sys.stderr, "Sanity check done. Resuming ICE job."
    icec.run()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("pickle_filename", help="Last successful pickle (ex: clusterOut/output/input.split_00.fa.pickle)")
    parser.add_argument("--root_dir", default="clusterOut/", help="root dir (default: clusterOut/)")
    parser.add_argument("--ccs_fofn", default="reads_of_insert.fofn", help="(default: reads_of_insert.fofn)")
    parser.add_argument("--flnc", default="isoseq_flnc.fasta", help="(default: isoseq_flnc.fasta)")
    parser.add_argument("--fasta_files_to_add", default=None, help="Comma-separated additional fasta files to add (default: None)")

    args = parser.parse_args()

    pickup_icec_job(args.pickle_filename, args.ccs_fofn, args.flnc, args.fasta_files_to_add, args.root_dir)
