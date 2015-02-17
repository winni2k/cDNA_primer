__author__ = 'etseng@pacificbiosciences.com'
import os, sys
from collections import defaultdict
from cPickle import dump, load
from csv import DictReader
from pbcore.io.FastaIO import FastaReader

def get_roi_len(seqid):
    # <movie>/<zmw>/<start>_<end>_CCS
    if not seqid.endswith('_CCS'):
        print >> sys.stderr, "Sequence ID format must be <movie>/<zmw>/<start>_<end>_CCS! Abort!"
        sys.exit(-1)
    s, e, junk = seqid.split('/')[2].split('_')
    return abs(int(s)-int(e))

def read_group_filename(group_filename, is_cid=True, sample_prefixes=None):
    """
    Make the connection between partitioned results and final (ex: PB.1.1)
    The partitioned results could either be ICE cluster (ex: i1_c123) or RoIs

    Return: dict of seq_or_ice_cluster --> collapsed cluster ID
    """
    cid_info = {} # ex: i1 --> c123 --> PB.1.1, or None --> c123 --> PB.1.1
    if sample_prefixes is not None:
        for sample_prefix in sample_prefixes:
            cid_info[sample_prefix] = {}
    else:
        cid_info[None] = {}
    for line in open(group_filename):
        pbid, members = line.strip().split('\t')
        for cid in members.split(','):
            # ex: x is 'i1_c123/f3p0/123 or
            # m131116_014707_42141_c100591062550000001823103405221462_s1_p0/93278/31_1189_CCS
            if sample_prefixes is None:
                if is_cid: cid = cid.split('/')[0]
                cid_info[None][cid] = pbid
            else:
                if any(cid.startswith(sample_prefix + '|') for sample_prefix in sample_prefixes):
                    sample_prefix, cid = cid.split('|', 1)
                    if is_cid: cid = cid.split('/')[0]
                    cid_info[sample_prefix][cid] = pbid
    return cid_info

def output_read_count_IsoSeq_csv(cid_info, csv_filename, output_filename, output_mode='w'):
    """
    Given an Iso-Seq csv output file w/ format:

    cluster_id,read_id,read_type

    Turn into read_stats.txt format:

    id \t length \t is_fl \t stat \t pbid
    """
    mapped = {}  # nFL seq -> list of (sample_prefix, cluster) it belongs to

    if output_mode == 'w':
        f = open(output_filename, 'w')
        f.write("id\tlength\tis_fl\tstat\tpbid\n")
    elif output_mode == 'a':
        f = open(output_filename, 'a')
    else:
        raise Exception, "Output mode {0} not valid!".format(output_mode)

    for r in DictReader(open(csv_filename), delimiter=','):
        cid = str(r['cluster_id'])
        assert cid.startswith('c')
        x = r['read_id']
        if cid in cid_info:
            # if is FL, must be unique
            if r['read_type'] == 'FL':
                pbid, stat = cid_info[cid], 'unique'
                f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                    id=x, len=get_roi_len(x), is_fl='Y', stat=stat, pbid=pbid))
            else: # nonFL could be multi-mapped, must wait and see
                assert r['read_type'] == 'NonFL'
                # is only potentially unmapped, add all (movie-restricted) members to unmapped holder
                pbid = cid_info[cid]
                if x not in mapped: mapped[x] = set()
                mapped[x].add(pbid)
        else:
            # unmapped
            pbid, stat = 'NA', 'unmapped'
            f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                id=x, len=get_roi_len(x), is_fl='Y', stat=stat, pbid=pbid))

    # now we can go through the list of mapped to see which are uniquely mapped which are not
    for seqid, pbids in mapped.iteritems():
        if len(pbids) == 1: # unique
            stat = 'unique'
        else:
            stat = 'ambiguous'
        for pbid in pbids:
            f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                id=seqid, len=get_roi_len(seqid), is_fl='N', stat=stat, pbid=pbid))

    f.close()


def output_read_count_RoI(cid_info, roi_filename, output_filename):
    """
    For each
    """
    f = open(output_filename, 'w')
    f.write("id\tlength\tis_fl\tstat\tpbid\n")
    for r in FastaReader(roi_filename):
        if r.id in cid_info: pbid, stat = cid_info[r.name], 'unique'
        else: pbid, stat = 'NA', 'unmapped'
        f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
            id=r.name, len=get_roi_len(r.name), is_fl='Y', stat=stat, pbid=pbid))
    f.close()

def output_read_count_FL(cid_info, pickle_prefix_list, output_filename, output_mode='w', restricted_movies=None):
    """
    If restricted_movies is None, all nonFL reads are output.
    Otherwise (esp. in case where I binned by size then pooled in the end), give the list of movies associated
    with a particular list of cell runs (ex: brain_2to3k_phusion nonFL only)

    Because may have multiple pickles, can ONLY determine which FL reads are unmapped at the VERY END
    """
    unmapped_holder = set() # will hold anything that was unmapped in one of the pickles
    mapped_holder = set() # will hold anything that was mapped in (must be exactly) one of the pickles
    # then to get the true unmapped just to {unmapped} - {mapped}

    if output_mode == 'w':
        f = open(output_filename, 'w')
        f.write("id\tlength\tis_fl\tstat\tpbid\n")
    elif output_mode == 'a':
        f = open(output_filename, 'a')
    else:
        raise Exception, "Output mode {0} not valid!".format(output_mode)

    for sample_prefix, pickle_filename in pickle_prefix_list:
        with open(pickle_filename) as h:
            uc = load(h)['uc']
        for cid_no_prefix, members in uc.iteritems():
            cid = 'c' + str(cid_no_prefix)
            if cid in cid_info[sample_prefix]:
                # can immediately add all (movie-restricted) members to mapped
                pbid, stat = cid_info[sample_prefix][cid], 'unique'
                for x in members:
                    movie = x.split('/')[0]
                    if restricted_movies is None or movie in restricted_movies:
                        mapped_holder.add(x)
                        f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                        id=x, len=get_roi_len(x), is_fl='Y', stat=stat, pbid=pbid))
            else:
                # is only potentially unmapped, add all (movie-restricted) members to unmapped holder
                for x in members:
                    movie = x.split('/')[0]
                    if restricted_movies is None or movie in restricted_movies:
                        unmapped_holder.add(x)

    # now with all the pickles processed we can determine which of all (movie-restricted) FL reads
    # are not mapped in any of the pickles
    unmapped_holder = unmapped_holder.difference(mapped_holder)
    pbid, stat = 'NA', 'unmapped'
    for x in unmapped_holder:
        f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
        id=x, len=get_roi_len(x), is_fl='Y', stat=stat, pbid=pbid))
    f.close()

def output_read_count_nFL(cid_info, pickle_prefix_list, output_filename, output_mode='w', restricted_movies=None):
    """
    If restricted_movies is None, all nonFL reads are output.
    Otherwise (esp. in case where I binned by size then pooled in the end), give the list of movies associated
    with a particular list of cell runs (ex: brain_2to3k_phusion nonFL only)

    There is no guarantee that the non-FL reads are shared between the pickles, they might be or not
    Instead determine unmapped (movie-restricted) non-FL reads at the very end
    """
    unmapped_holder = set() # will hold anything that was unmapped in one of the pickles
    # then to get the true unmapped just to {unmapped} - {mapped}

    mapped = {}  # nFL seq -> list of (sample_prefix, cluster) it belongs to

    if output_mode == 'w':
        f = open(output_filename, 'w')
        f.write("id\tlength\tis_fl\tstat\tpbid\n")
    elif output_mode == 'a':
        f = open(output_filename, 'a')
    else:
        raise Exception, "Output mode {0} not valid!".format(output_mode)

    for sample_prefix, pickle_filename in pickle_prefix_list:
        with open(pickle_filename) as h:
            result = load(h)
            uc = result['partial_uc']
            if restricted_movies is None:
                unmapped_holder.update(result['nohit'])
            else:
                unmapped_holder.update(filter(lambda x: x.split('/')[0] in restricted_movies, result['nohit']))

        for cid_no_prefix, members in uc.iteritems():
            cid = 'c' + str(cid_no_prefix)
            if cid in cid_info[sample_prefix]: # is at least mapped
                pbid = cid_info[sample_prefix][cid]
                for x in members:
                    movie = x.split('/')[0]
                    if restricted_movies is None or movie in restricted_movies:
                        if x not in mapped: mapped[x] = set()
                        mapped[x].add(pbid)
            else: # not entirely sure it is unmapped but put it in in the meantime
                for x in members:
                    movie = x.split('/')[0]
                    if restricted_movies is None or movie in restricted_movies:
                        unmapped_holder.add(x)

    # now we can go through the list of mapped to see which are uniquely mapped which are not
    for seqid, pbids in mapped.iteritems():
        if len(pbids) == 1: # unique
            stat = 'unique'
        else:
            stat = 'ambiguous'
        for pbid in pbids:
            f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                id=seqid, len=get_roi_len(seqid), is_fl='N', stat=stat, pbid=pbid))

    unmapped_holder = unmapped_holder.difference(mapped)

    # write the nohits
    for seqid in unmapped_holder:
        f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
            id=seqid, len=get_roi_len(seqid), is_fl='N', stat='unmapped', pbid='NA'))
    f.close()


def make_abundance_file(read_count_filename, output_filename, given_total=None, restricted_movies=None, write_header_comments=True):
    """
    If given_total is not None, use it instead of the total count based on <read_count_filename>
    given_total should be dict of {fl, nfl, nfl_amb}
    """
    total_ids = {'fl':set(), 'nfl':set(), 'nfl_amb':set()}
    tally = defaultdict(lambda: {'fl':0, 'nfl':0, 'nfl_amb':0}) # pbid, could be NA --> # of FL reads mapped to it
    amb_count = defaultdict(lambda: []) # non-fl id --> list of pbid matches

    reader = DictReader(open(read_count_filename), delimiter='\t')
    for r in reader:
        movie = r['id'].split('/')[0]
        if restricted_movies is None or movie in restricted_movies:
            if r['pbid'] != 'NA':
                if r['is_fl'] == 'Y': # FL, must be uniquely mapped
                    assert r['stat'] == 'unique'
                    tally[r['pbid']]['fl'] += 1
                    total_ids['fl'].add(r['id'])
                else: # non-FL, can be ambiguously mapped
                    if r['stat'] == 'unique':
                        tally[r['pbid']]['nfl'] += 1
                        total_ids['nfl'].add(r['id'])
                    else:
                        assert r['stat'] == 'ambiguous'
                        amb_count[r['id']].append(r['pbid'])
                        total_ids['nfl_amb'].add(r['id'])
            else: # even if it is unmapped it still counts in the abundance total!
                if r['is_fl'] == 'Y':
                    total_ids['fl'].add(r['id'])
                else:
                    total_ids['nfl'].add(r['id'])


    # put the ambiguous back in tally weighted
    for seqid, pbids in amb_count.iteritems():
        weight = 1. / len(pbids)
        for pbid in pbids:
             tally[pbid]['nfl_amb'] += weight

    if given_total is not None:
        use_total_fl = given_total['fl']
        use_total_nfl = given_total['fl'] + given_total['nfl']
        # ToDo: the below is NOT EXACTLY CORRECT!! Fix later!
        use_total_nfl_amb = given_total['fl'] + given_total['nfl'] + given_total['nfl_amb']
    else:
        use_total_fl = len(total_ids['fl'])
        use_total_nfl = len(total_ids['fl']) + len(total_ids['nfl'])
        use_total_nfl_amb = len(total_ids['fl']) + len(total_ids['nfl']) + len(total_ids['nfl_amb'])

    f = open(output_filename,'w')
    if write_header_comments:
        f.write("#\n")
        f.write("# -----------------\n")
        f.write("# Field explanation\n")
        f.write("# -----------------\n")
        f.write("# count_fl: Number of associated FL reads\n")
        f.write("# count_nfl: Number of associated FL + unique nFL reads\n")
        f.write("# count_nfl_amb: Number of associated FL + unique nFL + weighted ambiguous nFL reads\n")
        f.write("# norm_fl: count_fl / total number of FL reads\n")
        f.write("# norm_nfl: count_nfl / total number of FL + unique nFL reads\n")
        f.write("# norm_nfl_amb: count_nfl_amb / total number of all reads\n")
        f.write("# Total Number of FL reads: {0}\n".format(use_total_fl))
        f.write("# Total Number of FL + unique nFL reads: {0}\n".format(use_total_nfl))
        f.write("# Total Number of all reads: {0}\n".format(use_total_nfl_amb))
        f.write("#\n")
    f.write("pbid\tcount_fl\tcount_nfl\tcount_nfl_amb\tnorm_fl\tnorm_nfl\tnorm_nfl_amb\n")

    keys = tally.keys()
    keys.sort(key=lambda x: map(int, x.split('.')[1:])) # sort by PB.1, PB.2....
    for k in keys:
        v = tally[k]
        a = v['fl']
        b = a + v['nfl']
        c = b + v['nfl_amb']
        f.write("{0}\t{1}\t{2}\t{3}\t".format(k, a, b, c))
        f.write("{0:.4e}\t{1:.4e}\t{2:.4e}\n".format(a*1./use_total_fl, b*1./use_total_nfl, c*1./use_total_nfl_amb))
    f.close()


