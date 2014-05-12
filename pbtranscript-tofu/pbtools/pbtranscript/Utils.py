"""Define util functions."""
import os.path as op
import os
import shutil
import logging
import sys

def revcmp(seq):
    """Given a sequence return its reverse complement sequence."""
    NTMAP = {'a':'t', 'c':'g', 't':'a', 'g':'c',
             'A':'T', 'C':'G', 'T':'A', 'G':'C'}
    return "".join([NTMAP[x] for x in seq])[::-1]


def realpath(f):
    """Return absolute, user expanded path."""
    if f is None:
        return None
    return op.abspath(op.expanduser(f))


def mkdir(path):
    """Create a directory if it does not pre-exist,
    otherwise, pass."""
    if not op.exists(path):
        os.makedirs(path)


def mknewdir(path):
    """Create a new directory if it does not pre-exist,
    otherwise, delete it and then re-create it."""
    if op.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def touch(path):
    """touch a file."""
    if op.exists(path):
        os.utime(path, None)
    else:
        open(path, 'a').close()


def generateChunkedFN(out_dir, prefix, num_chunks):
    """Generate n chunked file names, e.g.
    outDir/$prefix.0, outDir/$prefix.1, ..., outDir/$prefix.num_chunks-1
    """
    return [op.join(out_dir, prefix + "." + str(i))
            for i in xrange (0, num_chunks)]


def get_files_from_fofn(fofn_filename):
    """Return a list of file names within a fofn file."""
    fns = []
    try:
        with open(fofn_filename, 'r') as fofn:
            for line in fofn:
                fns.append(realpath(line.strip()))
    except (IOError, OSError) as e:
        raise IOError("Failed to read from fofn file {fofn}.\n".
                      format(fofn=fofn_filename) + str(e))
    return fns


def write_files_to_fofn(file_names, fofn_filename):
    """Write files in list file_names to fofn_filename."""
    try:
        with open(fofn_filename, 'w') as fofn:
            for fn in file_names:
                fofn.write(str(fn) + "\n")
    except (IOError, OSError) as e:
        raise IOError("Failed to files to fofn file {fofn}.\n".
                      format(fofn=fofn_filename) + str(e))


def setup_log(alog, file_name=None, level=logging.DEBUG, str_formatter=None):
    """
    Copied from mkocher's pbreports/utils.py.
    Util function for setting up logging.

    Due to how smrtpipe logs, the default behavior is that the stdout
    is where the logging is redirected. If a file name is given the log
    will be written to that file.

    :param log: (log instance) Log instance that handlers and filters will
    be added.
    :param file_name: (str, None), Path to file. If None, stdout will be used.
    :param level: (int) logging level
    """
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)

    if str_formatter is None:
        str_formatter = '[%(levelname)s] %(asctime)-15s ' + \
                        '[%(name)s %(funcName)s %(lineno)d] %(message)s'

    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    alog.addHandler(handler)
    alog.setLevel(level)


def now_str():
    """Return string of current time."""
    import datetime
    return str(datetime.datetime.now()).split(".")[0]


def phred_to_qv(phred):
    """Phred value to quality value."""
    return 10 ** -(phred/10.0)


def cat_files(src, dst):
    """Concatenate files in src and save to dst.
       src --- source file names in a list
       dst --- destinate file name
    """
    if src is None or len(src) == 0:
        raise ValueError("src should contain at least one file.")
    if dst in src:
        raise IOError("Unable to cat a file and save to itself.")

    with open (dst, 'w') as writer:
        for src_f in src:
            with open(src_f, 'r') as reader:
                for line in reader:
                    writer.write(line.rstrip() + '\n')


def get_all_files_in_dir(dir_path, extension=None):
    """return all files in a directory."""
    fs = []
    for f in os.listdir(dir_path):
        if extension is None:
            fs.append(f)
        else:
            if f.endswith(extension):
                fs.append(f)
    return fs


class CIGAR(object):
    """Cigar string."""
    def __init__(self, cigar_str):
        self.cigar_str = cigar_str
        self.num_match = 0
        self.num_mismatch = 0
        self.num_insert = 0
        self.num_deletion = 0
        self.num_hardclip = 0
        self.num_softclip = 0
        self.num_padding = 0
        self.num_unknown = 0
        self.parse(self.cigar_str)

    def parse(self, cigar_str):
        """Parse cigar string."""
        s = cigar_str
        if s == "*" or s == "=":
            return
        while(len(s) > 0):
            i = 0
            while s[i].isdigit() and i < len(s):
                i += 1
            num = int(s[:i])
            action = s[i]
            s = s[i+1:]
            if action == "M":
                self.num_match += num
            elif action == "X":
                self.num_mismatch += num
            elif action == "I":
                self.num_insert += num
            elif action == "D":
                self.num_deletion += num
            elif action == "H":
                self.num_hardclip += num
            elif action == "S":
                self.num_softclip += num
            elif action == "P":
                self.num_padding += num
            elif action == "N":
                self.num_unknown += num
            else:
                raise ValueError("Can not parse CIGAR string " +
                                 "{s}".format(s=cigar_str))
# using regular expression is 20% slower than naive way
#        pattern = r"^(\d+)(M|I|D|X|S|H|P|N)(.*)$"
#        while(len(s) > 0):
#            m = re.search(pattern, s)
#            if m:
#                num = int(m.groups()[0])
#                action = m.groups()[1]

    def __str__(self):
        return "Match = {m}, ".format(m=self.num_match) + \
               "Mismatch = {m}, ".format(m=self.num_mismatch) + \
               "Insert = {m}, ".format(m=self.num_insert) + \
               "Deletion = {m}, ".format(m=self.num_deletion) + \
               "HardClipping = {m}, ".format(m=self.num_hardclip) + \
               "SoftClipping = {m}, ".format(m=self.num_softclip) + \
               "Padding = {m}, ".format(m=self.num_padding)

    def match_seq(self, seq):
        """Return if this cigar string matches the sequence."""
        return self.cigar_str == "*" or self.cigar_str == "=" or \
            (self.num_match + self.num_insert + self.num_softclip == len(seq))


def cigar_match_seq(sam_str):
    """Return True if cigar length match sequence length, otherwise, False"""
    fields = sam_str.split('\t')
    cigar, seq = CIGAR(fields[5]), fields[9]
    if cigar.match_seq(seq):
        return True
    else:
        return False


def filter_sam(in_sam, out_sam):
    """Filter sam alignments with bad cigar string."""
    if not op.exists(in_sam):
        raise IOError("Unable to find input sam {f}".format(f=in_sam))
    if realpath(in_sam) == realpath(out_sam):
        raise IOError("in_sam and out_sam can not be identical.")
    with open(in_sam, 'r') as reader, \
         open(out_sam, 'w') as writer:
        for l, in_line in enumerate(reader):
            logging.info("Processing {l}".format(l=in_line))
            if in_line.startswith("#") or \
               in_line.startswith("@") or \
               cigar_match_seq(in_line):
                writer.write(in_line)
            else:
                logging.warn("line {l}, cigar does not match sequence.".
                             format(l=l+1))


def ln(src, dst):
    """if src and dst are identical, pass. Otherwise, create dst, a soft
    symbolic link pointing to src."""
    if realpath(src) !=  realpath(dst):
        if op.exists(dst):
            os.remove(dst)
        logging.debug("Creating a symbolic link {dst} pointing to {src}".
                     format(dst=dst, src=src))
        os.symlink(src, dst)

