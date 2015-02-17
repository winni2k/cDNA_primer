#!/usr/bin/env python
"""
Find initial mutual exclusive cliques by aligning
input reads against itself.
"""
__author__ = 'etseng@pacificbiosciences.com'

import os
import networkx as nx
import logging
from pbcore.io import FastaReader
from pbcore.util.Process import backticks
import pbtools.pbtranscript.ice.pClique as pClique
from pbtools.pbtranscript.ice.IceUtils import blasr_against_ref


class IceInit(object):
    """Iterative clustering and error correction."""
    def __init__(self, readsFa, qver_get_func,
        ice_opts, sge_opts):

        self.readsFa = readsFa
        self.qver_get_func = qver_get_func
        self.ice_opts = ice_opts
        self.sge_opts = sge_opts
        self.uc = None

        self.uc = self.init_cluster_by_clique(
            readsFa=readsFa,
            qver_get_func=qver_get_func,
            ice_opts=ice_opts, sge_opts=sge_opts)

    def _align(self, queryFa, targetFa, outFN, ice_opts, sge_opts):
        """Align input reads against itself using BLASR."""
        if os.path.exists(outFN):
            logging.info("{0} already exists. No need to run BLASR.".
                format(outFN))
        else:
            cmd = "blasr {q} {t} ".format(q=queryFa, t=targetFa) + \
                  "-m 5 -maxLCPLength 15 " + \
                  "-nproc {cpu} ".format(cpu=sge_opts.blasr_nproc) + \
                  "-maxScore {score} ".format(score=ice_opts.maxScore) + \
                  "-bestn {n} -nCandidates {n} ".format(n=ice_opts.bestn) + \
                  "-out {o}".format(o=outFN)
            logging.info("Calling {cmd}".format(cmd=cmd))
            _output, code, msg = backticks(cmd)
            if code != 0:
                errMsg = "{cmd} exited with {code}: {msg}".\
                        format(cmd=cmd, code=code, msg=msg)
                logging.error(errMsg)
                raise RuntimeError (errMsg)

    def _makeGraphFromM5(self, m5FN, qver_get_func, ice_opts):
        """Construct a graph from a BLASR M5 file."""
        alignGraph = nx.Graph()

        for r in blasr_against_ref(output_filename=m5FN,
            is_FL=True,
            sID_starts_with_c=False,
            qver_get_func=qver_get_func,
            ece_penalty=ice_opts.ece_penalty,
            ece_min_len=ice_opts.ece_min_len):
            if r.qID == r.cID:
                continue # self hit, ignore
            if r.ece_arr is not None:
                logging.debug("adding edge {0},{1}".format(r.qID, r.cID))
                alignGraph.add_edge(r.qID, r.cID)
        return alignGraph

    def _findCliques(self, alignGraph, readsFa):
        """
        Find all mutually exclusive cliques within the graph, with decreased
        size.

        alignGraph - a graph, each node represent a read and each edge
        represents an alignment between two end points.

        Return a dictionary of clique indices and nodes.
            key = index of a clique
            value = nodes within a clique
        Cliques are ordered by their size descendingly: index up, size down
        Reads which are not included in any cliques will be added as cliques
        of size 1.
        """
        uc = {}    # To keep cliques found
        used = []  # nodes within any cliques
        ind = 0    # index of clique to discover

        deg = alignGraph.degree().items()
        # Sort tuples of (node, degree) by degree, descendingly
        deg.sort(key=lambda x:x[1], reverse=True)
        for d in deg:
            node = d[0]  # node which has the largest degree in alignGraph
            if node not in alignGraph:
                continue
            # just get the immediate neighbors since we're looking for perfect
            # cliques
            subGraph = alignGraph.subgraph([node] + alignGraph.neighbors(node))
            subNodes = subGraph.nodes()
            # Convert from networkx.Graph to a sparse matrix
            S, H = pClique.convert_graph_connectivity_to_sparse(
                subGraph, subNodes)
            # index of the 'node' in the sub-graph
            seed_i = subNodes.index(node)
            # Grasp a clique from subGraph, and return indices of clique nodes
            tQ = pClique.grasp(S, H, 1., 5, seed_i)
            if len(tQ) > 0:
                c = [subNodes[i] for i in tQ]  # nodes in the clique
                uc[ind] = c  # Add the clique to uc
                ind += 1
                used += c    # Add clique nodes to used
                # Remove clique nodes from alignGraph and continue
                alignGraph.remove_nodes_from(c)

        with FastaReader(readsFa) as reader:
            for r in reader:
                rid = r.name.split()[0]
                if rid not in used:
                    uc[ind] = [rid]
                    ind += 1
        return uc

    def init_cluster_by_clique(self, readsFa, qver_get_func,
            ice_opts, sge_opts):
        """
        Only called once and in the very beginning, when (probably a subset)
        of sequences are given to generate the initial cluster.

        readsFa --- initial fasta filename, probably called *_split00.fa
        qver_get_func --- function that returns QVs on reads
        bestn --- parameter in BLASR, higher helps in finding perfect
            cliques but bigger output
        nproc, maxScore --- parameter in BLASR, set maxScore appropriate
            to input transcript length
        ece_penalty, ece_min_len --- parameter in isoform hit calling

        Self-blasr input then iteratively find all mutually exclusive
            cliques (in decreasing size)
        Returns dict of cluster_index --> list of seqids
        which is the 'uc' dict that can be used by IceIterative
        """
        outFN = readsFa + '.self.blasr'

        self._align(queryFa=readsFa, targetFa=readsFa, outFN=outFN,
            ice_opts=ice_opts, sge_opts=sge_opts)

        alignGraph = self._makeGraphFromM5(m5FN=outFN,
            qver_get_func=qver_get_func, ice_opts=ice_opts)

        uc = self._findCliques(alignGraph=alignGraph, readsFa=readsFa)
        return uc

