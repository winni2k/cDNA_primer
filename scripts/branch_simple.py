import os, sys
import pdb
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import BioReaders
import c_branch
from intersection_unique import IntervalTreeUnique

INTINF = 9999


class ContiVec:
    """
    Original struct: 'BC'

    """
    def __init__(self, size):
        self.baseC = np.zeros(size, dtype=np.int) # this was .B in the original code, base coverage
        self.altC = np.zeros(size, dtype=np.int)  # this was .C in the original code, evidence for alternative junction
        self.transC = np.zeros(size, dtype=np.int)  # this is .isTrans in the original code


class BranchSimple:
    """
    BranchSimple is designed for just creating exons from PacBio's GMAP results
    Does not use Illumina
    The genome fasta file is only used to figure out the size of genomes
    """
    def __init__(self, transfrag_filename, cov_threshold=2, cuff_index=1):
        self.contiVec = None # current ContiVec object
        self.exons = None

        self.transfrag_filename = transfrag_filename
        self.transfrag_len_dict = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(transfrag_filename), 'fasta'))

        self.cov_threshold = cov_threshold # only output GTF records if >= this many GMAP records support it (this must be if I'm running non-clustered fasta on GMAP)

        self.cuff_index = cuff_index

    def run_on_gmap_sam(self, gmap_sam_filename, output_prefix):
        f_good = open(output_prefix + '.transcripts.gtf', 'w')
        f_bad = open(output_prefix + '.skipped.gtf', 'w')
        for recs in self.iter_gmap_sam(gmap_sam_filename):
            if len(recs['+']) > 0: self.process_records(recs['+'], f_good, f_bad)
            if len(recs['-']) > 0: self.process_records(recs['-'], f_good, f_bad)
        f_good.close()
        f_bad.close()


    def iter_gmap_sam(self, gmap_sam_filename):
        """
        Iterate over a SORTED GMAP SAM file.
        Return a collection of records that overlap by at least 1 base.
        """
        def sep_by_strand(records):
            output = {'+':[], '-':[]}
            for r in records:
                output[r.flag.strand].append(r)
            return output

        records = None # holds the current set of records that overlap in coordinates
        iter = BioReaders.GMAPSAMReader(gmap_sam_filename, True, query_len_dict=self.transfrag_len_dict)
        for r in iter:
            if r.sID != '*': break
        records = [r]
        max_end = r.sEnd

        for r in iter:
            if r.sID == records[0].sID and r.sStart < records[-1].sStart:
                print >> sys.stderr, "SAM file is NOT sorted. ABORT!"
                sys.exit(-1)
            if r.sID != records[0].sID or r.sStart > max_end:
                yield sep_by_strand(records)
                records = [r]
                max_end = r.sEnd
            else:
                records.append(r)
                max_end = max(max_end, r.sEnd)
        yield sep_by_strand(records)


    def parse_transfrag2contig(self, gmap_sam_records):
        """
        GMAP SAM file MUST BE SORTED! (same criterion as cufflinks)

        Goes through
        """
        records = gmap_sam_records
        # first figure out how long the "pseudo-chromosome" size is
        offset = records[0].sStart
        self.offset = offset
        self.chrom = records[0].sID
        self.strand = records[0].flag.strand
        chrom_size = max(x.sEnd for x in records) - records[0].sStart
        self.contiVec = ContiVec(chrom_size)
        for r in records:
            for i,e in enumerate(r.segments):
                # fill base coverage
                self.contiVec.baseC[(e.start-offset):(e.end-offset)] += 1

                # in the original code, the mapped start altC was set to -MAX and end to MAX
                if i != 0: # skip the alt if this is beginning of first exon
                    self.contiVec.altC[(e.start-offset)] -= INTINF
                if i != len(r.segments)-1: # skip if end of last exon
                    self.contiVec.altC[(e.end-offset-1)] += INTINF  # adjust to 0-based
                # set .transC
                self.contiVec.transC[(e.start-offset):(e.end-offset)] += 1


    def exon_finding(self):
        """
        Go through contiVec to identify the exons using base coverage (.baseC) and alt junction evidence (.altC)
        """
        v = self.contiVec
        self.exons = c_branch.exon_finding(v.baseC, v.altC, v.transC, \
                            len(v.baseC), 10, 0, self.offset)


    def match_record(self, r, tolerate_middle=10, tolerate_end=10**4, ok_to_miss_matches=False, intervals_adjacent=True):
        """
        r --- a gmap sam record
        """
        result = []
        num_exons = len(r.segments)
        for i,e in enumerate(r.segments):
            # allow the first and last exon to be longer or shorter in either the long read or the detected exons
            # however for middle exons the matching should be pretty precise
            # ToDo: parameterize this
            if i == 0 and num_exons > 1:
                tolerate_l = tolerate_end
                tolerate_r = tolerate_middle
            elif i == len(r.segments)-1 and num_exons > 1:
                tolerate_l = tolerate_middle
                tolerate_r = tolerate_end
            else:
                tolerate_l = tolerate_middle
                tolerate_r = tolerate_middle
            matches = c_branch.exon_matching(self.exons, e, tolerate_l, tolerate_r, intervals_adjacent) # ToDo: CHANGE this back to c_branch later
            if matches is None:
                if not ok_to_miss_matches:
                    return None
            else:
                if (len(result) >= 1 and result[-1].value >= matches[0].value) and (not ok_to_miss_matches):
                    return None
                result += matches
        return result


    def process_records(self, records, f_good, f_bad):
        """
        Given a set of records
        (1) process them by running through parse_transfrag2contig
        (2) call exons by exon_finding
        (3) go through each record, get the list of "nodes" they corresspond to
        (4) collapse identical records (must have same set of node values)

        Write out to GTF format
        """
        self.parse_transfrag2contig(records)
        self.exon_finding()

        bad_index = -1

        cluster = defaultdict(lambda: 0)
        for r in records:
            stuff = self.match_record(r)
            if stuff is None: # ToDo: handle this better later
                # write this out to f_bad
                f_bad.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"PB.{i}.{j}\"; cov \"{cov}\";\n".format(\
                    chr=self.chrom, s=r.segments[0].start+1, e=r.segments[-1].end, i=self.cuff_index, j=bad_index, strand=self.strand, cov=1))
                for seg in r.segments:
                    f_bad.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"PB.{i}.{j}\"; cov \"{cov}\";\n".format(\
                    chr=self.chrom, s=seg.start+1, e=seg.end, i=self.cuff_index, j=bad_index, strand=self.strand, cov=1))
                bad_index -= 1
            else:
                cluster[frozenset(x.value for x in stuff)] += 1

        self.isoform_index = 1
        # make the exon value --> interval dictionary
        a = []
        self.exons.traverse(a.append)
        node_d = {}
        for x in a: node_d[x.interval.value] = x

        keys = list(cluster.keys())
        keys.sort()
        for node_set in keys:
            if cluster[node_set] < self.cov_threshold:
                f_out = f_bad
                index = bad_index
                bad_index -= 1
            else:
                f_out = f_good
                index = self.isoform_index
                self.isoform_index += 1
            node_set_list = list(node_set)
            node_set_list.sort()
            segments = [node_d[x] for x in node_set_list]
            f_out.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"PB.{i}.{j}\"; cov \"{cov}\";\n".format(\
                chr=self.chrom, s=segments[0].start+1, e=segments[-1].end, i=self.cuff_index, j=index, strand=self.strand, cov=cluster[node_set]))

            i = 0
            j = 0
            for j in xrange(1, len(segments)):
                if segments[j].start != segments[j-1].end:
                    f_out.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"PB.{i}.{j}\"; cov \"{cov}\";\n".format(\
                    chr=self.chrom, s=segments[i].start+1, e=segments[j-1].end, i=self.cuff_index, j=index, strand=self.strand, cov=cluster[node_set]))
                    i = j
            f_out.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"PB.{i}.{j}\"; cov \"{cov}\";\n".format(\
                    chr=self.chrom, s=segments[i].start+1, e=segments[j].end, i=self.cuff_index, j=index, strand=self.strand, cov=cluster[node_set]))

        self.cuff_index += 1

        return cluster


