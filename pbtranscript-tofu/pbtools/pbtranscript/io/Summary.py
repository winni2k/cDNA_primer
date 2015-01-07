"""Define ClassifySummary and ClusterSummary."""


class Summary(object):

    """Supert class for ClassifySummary and ClusterSummary."""
    @property
    def fieldsNames(self):
        """Specific for Classify and Cluster."""
        raise NotImplementedError("Summary.fieldsNames() not implemented")

    @property
    def fields(self):
        """Specific for Classify and Cluster."""
        raise NotImplementedError("Summary.fields() not implemented")

    def __str__(self):
        return "\n".join(["{name}={val}".format(name=name, val=val)
                          for name, val in zip(self.fieldsNames, self.fields)])

    def write(self, outFile):
        """Write summary to outFile."""
        assert(len(self.fieldsNames) == len(self.fields))
        with open(outFile, 'w') as writer:
            writer.write(self.__str__() + "\n")


class ClassifySummary(Summary):

    """A ClassifySummary object has all classify summary attributes."""

    def __init__(self):
        Summary.__init__(self)

        self.num_reads = 0  # number of reads from in.fasta
        self.num_5_seen = 0  # number 5' primer seen reads within in.fasta
        self.num_3_seen = 0  # number 3' primer seen reads within in.fasta
        self.num_polyA_seen = 0      # number polyA seen reads within in.fasta

        # number of filtered short reads whose length < min_seq_len
        self.num_filtered_short_reads = 0

        # number of full-length reads, > min_seq_len,
        # can either be chimeric or non-chimeric
        self.num_fl = 0
        # number of full-length non-chimeric reads, > min_seq_len
        self.num_flnc = 0
        # number of full-length chimeric reads
        self.num_flc = 0
        # total number of bases in full-length non-chimeric reads
        self.num_flnc_bases = 0

        # number of non-full-length reads, > min_seq_len,
        # can either be chimeric or non-chimeric
        self.num_nfl = 0
        # number of non-full-length non-chimeric reads
        self.num_nflnc = None
        # number of non-full-length chimeric reads
        self.num_nflc = None

    @property
    def avg_flnc_len(self):
        """Return average read length of full-length non-chimeric reads."""
        return int(self.num_flnc_bases / self.num_flnc)

    @property
    def fieldsNames(self):
        """Return all fields names in a list.
        These fields will be displayed as attributes in json files, SMRTPipe
        and SMRTPortal. Changes to field names will lead to changes of
        //depot/software/assembly/java/smrtportal/src/main/resources/reportRules/*/
        RS_IsoSeq.1.rules.xml, isoseqClassifyRules.xml and isoseqClusterRules.xml.
        """
        rets = ["Number of reads of insert",
                "Number of five prime reads",
                "Number of three prime reads",
                "Number of poly-A reads",
                "Number of filtered short reads",
                "Number of non-full-length reads",
                "Number of full-length reads",
                "Number of full-length non-chimeric reads",
                "Average full-length non-chimeric read length"]
        if self.num_nflnc is not None and self.num_nflc is not None:
            rets.extend(["Number of non-full-length non-chimeric reads",
                         "Number of non-full-length chimeric reads"])
        return rets

    @property
    def fields(self):
        """Return fiels values in a list."""
        rets = [self.num_reads,
                self.num_5_seen,
                self.num_3_seen,
                self.num_polyA_seen,
                self.num_filtered_short_reads,
                self.num_nfl,
                self.num_fl,
                self.num_flnc,
                self.avg_flnc_len]
        if self.num_nflnc is not None and self.num_nflc is not None:
            rets.extend([self.num_nflnc, self.num_nflc])
        return rets


class ClusterSummary(Summary):

    """A ClusterSummary object has all cluster summary attributes."""

    def __init__(self):
        Summary.__init__(self)

        self.numConsensusIsoforms = 0  # number of consensus isoforms.
        # total number of bases in predicted consensus isoforms
        self.numTotalBases = 0
        self.num_polished_hq_isoforms = None
        self.num_polished_lq_isoforms = None

    @property
    def avgConsensusIsoformLength(self):
        """Return average read length of predicted consensus isoforms."""
        return int(self.numTotalBases / self.numConsensusIsoforms)

    @property
    def fieldsNames(self):
        """Return all fields names in a list.
        These fields will be displayed as attributes in json files, SMRTPipe
        and SMRTPortal. Changes to field names will lead to changes of
        //depot/software/assembly/java/smrtportal/src/main/resources/reportRules/*/
        RS_IsoSeq.1.rules.xml and isoseqClusterRules.xml.
        """
        fs = ["Number of consensus isoforms",
              "Average consensus isoforms read length"]
        if self.num_polished_hq_isoforms is not None:
            fs += ["Number of polished high-quality isoforms"]
        if self.num_polished_lq_isoforms is not None:
            fs += ["Number of polished low-quality isoforms"]
        return fs

    @property
    def fields(self):
        """Return fiels values in a list. Have to match self.fieldsNames"""
        fv = [self.numConsensusIsoforms,
              self.avgConsensusIsoformLength]
        if self.num_polished_hq_isoforms is not None:
            fv += [self.num_polished_hq_isoforms]
        if self.num_polished_lq_isoforms is not None:
            fv += [self.num_polished_lq_isoforms]
        return fv
