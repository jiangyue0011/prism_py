#!/usr/bin/env python

import sys
import pysam
import os
import commands
import argparse

class SV_cluster(object):
    def __init__(self, read):
        self.chr1 = read.reference_name
        self.chr2 = read.next_reference_name
        self.start = read.reference_start
        self.end = read.next_reference_start
        # the first read's mapping position in the cluster
        self.upstream_bnd = self.start
        # the last read's mapping position in the cluster
        self.downstream_bnd = self.end
        self.read_cnt = 1
        #self.length = abs(read.reference_start - read.next_reference_start) + read.query_length
        self.length = abs(read.reference_start - read.next_reference_start)
        if self.chr1 == self.chr2:
            if read.reference_start < read.next_reference_start:
                if not read.is_reverse and read.mate_is_reverse:
                    self.type = "DEL"
                elif read.is_reverse and not read.mate_is_reverse:
                    self.type = "DUP"
                else:
                    self.type = "INV"
            else:
                if not read.is_reverse and read.mate_is_reverse:
                    self.type = "DUP"
                elif read.is_reverse and not read.mate_is_reverse:
                    self.type = "DEL"
                else:
                    self.type = "INV"
        else:
            self.type = "CTX"

    def update(self, read):
        self.start = int(0.5*abs(self.start + read.reference_start))
        self.end = int(0.5*abs(self.end + read.next_reference_start))
        self.downstream_bnd = read.next_reference_start
        self.length = abs(self.start - self.end)
        self.read_cnt += 1

    def print_cluster(self):
        print self.chr1, self.start, self.chr2, self.end, self.type, self.read_cnt, self.upstream_bnd, self.downstream_bnd, self.length    

def pass_filter(read, mean, std):
    return True


dna_bam_file = sys.argv[1]
dna_bam =  pysam.AlignmentFile(dna_bam_file, "rb")
mean = 500
std = 30

cluster_del = []
cluster_inv = []
cluster_dup = []
cluster_ctx = []
clusters = {}
clusters.setdefault("DEL", cluster_del)
clusters.setdefault("INV", cluster_inv)
clusters.setdefault("DUP", cluster_dup)
clusters.setdefault("CTX", cluster_ctx)

read_cnt = 0
for read in dna_bam.fetch():
    if read.mapping_quality < 50:
        continue
    read_cnt += 1
    if read_cnt % 10000 == 0:
        print >> sys.stderr, read_cnt, "reads processed."
    read_cluster = SV_cluster(read)
    # skip normal pairs and reads whose mate have been processed
    if (read_cluster.length < mean + 5*std and read_cluster.type == "DEL") or (read_cluster.start > read_cluster.end):
        continue
    cluster_list = clusters[read_cluster.type]
    clustered = False
    bnd_idx = 0
    for idx,sv_cluster in enumerate(cluster_list):
        if read_cluster.chr1 == sv_cluster.chr1 and read_cluster.chr2 == sv_cluster.chr2 and abs(read_cluster.start - sv_cluster.start) < mean and abs(read_cluster.end - sv_cluster.end) < mean:
            sv_cluster.update(read)
            clustered = True
            break
        elif read_cluster.chr1 != sv_cluster.chr1 or read_cluster.start - sv_cluster.start > mean:
            bnd_idx = idx+1
        elif sv_cluster.start - read_cluster.start > mean:
            break
    if bnd_idx:
        for i in xrange(0, bnd_idx):
            if cluster_list[i].read_cnt > 2 and cluster_list[i].length > mean+5*std:
                cluster_list[i].print_cluster()
        del cluster_list[:bnd_idx]            
    if not clustered:
        cluster_list.append(read_cluster)

for sv_type in clusters:
    cluster_list = clusters[sv_type]
    for sv_cluster in cluster_list:
        if cluster_list[i].read_cnt > 2 and cluster_list[i].length > mean+5*std:
            sv_cluster.print_cluster()

    
    
    
