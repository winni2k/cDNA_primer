#!/usr/bin/env python
import os, sys, glob, shutil

def cleanup_ice(root_dir, delete_tmp_dir=False, delete_quiver_dir=False):
    """
    WARNING: should ONLY call this if SURE that ICE/Quiver is done!
    Will do some sanity checking to delete only if there's proof that ICE/Quiver is finished.
    """
    def del_files(files):
        for file in files:
            print >> sys.stderr, "Deleting file {0}....".format(file)
            os.remove(file)

    final_pickle = os.path.join(root_dir, 'output', 'final.pickle')
    if not os.path.exists(final_pickle):
        raise Exception, "{0} does not exist. ICE probably not finished!".format(final_pickle)

    files = glob.glob(os.path.join(root_dir, 'input.split*.fa*'))
    files += glob.glob(os.path.join(root_dir, 'ref_consensus*'))
    files += glob.glob(os.path.join(root_dir, 'tmp.orphan.*'))
    files += glob.glob(os.path.join(root_dir, 'current.fasta.*'))
    files += glob.glob(os.path.join(root_dir, '*dazz*'))
    files += glob.glob(os.path.join(root_dir, 'output', 'input.split*'))
    files += glob.glob(os.path.join(root_dir, 'output', 'tmp.consensus.*'))
    files = list(set(files))
    del_files(files)

    dirs_in_log = os.path.join(root_dir, 'log') # don't delete the log/*.log files
    # but delete things like log/0, log/1, log/2....
    for x in glob.glob(dirs_in_log):
        d = os.path.join(dirs_in_log, x)
        if os.path.isdir(d) and not d.endswith('.log'):
            shutil.rmtree(d)

    nfl_pickle = os.path.join(root_dir, 'output', 'map_noFL', 'nfl.all.partial_uc.pickle')
    if not os.path.exists(nfl_pickle):
        raise Exception, "{0} does not exist. Partial probably not finished!".format(nfl_pickle)

    files = glob.glob(os.path.join(root_dir, 'output', 'map_noFL', 'input.split*'))
    files += glob.glob(os.path.join(root_dir, 'output', 'map_noFL', '*dazz*'))
    files = list(set(files))
    del_files(files)

    if delete_tmp_dir:
        tmp_dir = os.path.join(root_dir, 'tmp')
        print >> sys.stderr, "Deleting {0}. Will take a while....".format(tmp_dir)
        shutil.rmtree(tmp_dir)

    if delete_quiver_dir:
        quiver_dir = os.path.join(root_dir, 'quivered')
        print >> sys.stderr, "Deleting {0}. Will take a while....".format(quiver_dir)
        shutil.rmtree(quiver_dir)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("root_dir")
    parser.add_argument("--delete_tmp", default=False, action="store_true", help="Delete cluster tmp dir (default: False)")
    parser.add_argument("--delete_quiver", default=False, action="store_true", help="Delete quivered dir (default: False)")

    args = parser.parse_args()

    cleanup_ice(args.root_dir, args.delete_tmp, args.delete_quiver)

