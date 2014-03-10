import re
from collections import namedtuple
from exceptions import StopIteration
       
class SimpleSAMReader:
    """
    A simplified SAM reader meant for speed. Skips CIGAR & FLAG parsing; identity/coverage calculation.
    """
    SAMheaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']    
    def __init__(self, filename, has_header):
        self.filename = filename
        self.f = open(filename)
        self.header = ''
        if has_header:
            while True:
                cur = self.f.tell()
                line = self.f.readline()
                if line[:3] not in SimpleSAMReader.SAMheaders:
                    break
                self.header += line
            self.f.seek(cur)
    
    def __iter__(self):
        return self
    
    def next(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return SimpleSAMRecord(line)    
    

def iter_cigar_string(cigar_string):
    num = cigar_string[0]
    for s in cigar_string[1:]:
        if str.isalpha(s):
            yield int(num), s
            num = ''
        else:
            num += s
    
class SimpleSAMRecord:
    cigar_rex = re.compile('(\d+)([MIDSHN])')
    SAMflag = namedtuple('SAMflag', ['is_paired', 'strand', 'PE_read_num'])
    def __init__(self, record_line):
        """
        Simple bare bones version: only has
        
        qID, sID, sStart, sEnd, qStart, qEnd, cigar
        
        Simplified assumptions:
        -- must be end-to-end alignment (so qStart always 0)
        -- must be unspliced (no 'N' in cigar string)
        """
        self.qID = None
        self.sID = None
        self.sStart = None
        self.sEnd = None
        self.qStart = 0
        self.qEnd = None # length of SEQ
        self.cigar = None

        self.process(record_line)

    def __str__(self):
        msg = \
        """
        qID: {q}
        sID: {s}
        sStart-sEnd: {ss}-{se}
        qStart-qEnd: {qs}-{qe}
        cigar: {c}
        """.format(q=self.qID, s=self.sID, \
            ss=self.sStart, se=self.sEnd, qs=self.qStart, qe=self.qEnd, c=self.cigar)
        return msg

    def parse_cigar(self, cigar, start):
        """
        M - match
        I - insertion w.r.t. to ref
        D - deletion w.r.t. to ref
        N - skipped (which means splice junction)
        S - soft clipped
        H - hard clipped (not shown in SEQ)

        ex: 50M43N3D

        NOTE: sets qStart & qEnd, which are often incorrect because of different ways to write CIGAR strings
              instead rely on XS/XE flags (from blasr or pbalign.py) to overwrite this later!!!

        Returns: genomic segment locations (using <start> as offset)
        """
        cur_end = start
        q_aln_len = 0
        #for num, type in SimpleSAMRecord.cigar_rex.findall(cigar):
        for num, type in iter_cigar_string(cigar):
            if type == 'I':
                q_aln_len += num
            elif type == 'M':
                cur_end += num
                q_aln_len += num
            elif type == 'D':
                cur_end += num
        self.qEnd = self.qStart + q_aln_len
        self.sEnd = cur_end

            
    def process(self, record_line):
        """
        """
        raw = record_line.split('\t')
        self.qID = raw[0]
        self.sID = raw[2]
        if self.sID == '*': # means no match! STOP here
            return
        self.sStart = int(raw[3]) - 1
        self.cigar = raw[5]
        self.parse_cigar(self.cigar, self.sStart)
        #self.flag = SimpleSAMRecord.parse_sam_flag(int(raw[1]))
        
        