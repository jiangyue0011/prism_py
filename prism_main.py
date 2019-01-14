#!/usr/bin/env python

import sys
import pysam

import smapper_core
import clusters
import utils

mean = 500
rlen = 100

def is_poorly_mapped(read):
    # soft-clipping or too many mismatch/indels
    if len(read.query_alignment_sequence) < read.query_length or int(read.get_tag("NM")) > 5:
        return True
    else:
        return False


def local2ref(local_start, ref1_start, ref1_len):
    if local_start > ref1_len:  # all bases aligned to ref2
        return -1
    else:
        return local_start + ref1_start


if __name__ == '__main__':
    ref_file = sys.argv[1]
    cluster_file = sys.argv[2]
    read_file = sys.argv[3]
    #read_file = sys.argv[2]
    # load all reference sequences into memory. Use indexed refrerence with pysam methods is slower.
    genome_seqs = smapper_core.Reference(ref_file)
    read_bam = pysam.AlignmentFile(read_file, "rb")
    cluster_list = clusters.load_clusters(cluster_file) # each element of cluster_list is a tuple, (ref, sv_start, sv_end, sv_type)
    for cluster in cluster_list:
        ref, sv_start, sv_end, sv_type = cluster
        #print sv_start-rlen, sv_start+0.5*mean
        ref1_start = sv_start - rlen
        ref1 = genome_seqs.ref_dict[ref][ref1_start:sv_start+mean/2]
        ref2 = genome_seqs.ref_dict[ref][sv_end-mean/2:sv_end+rlen]
        ref_overlap_len = min(0, (sv_end-mean)-(sv_start+mean))
        for read in read_bam.fetch(ref, sv_start-mean/2, sv_start+mean/2):
            #print read.query_name, is_poorly_mapped(read)
            if is_poorly_mapped(read):
                align_start, cigar, score =  smapper_core.sw_split(ref1, ref2, read.query_sequence, ref_overlap_len)
                ref_start = local2ref(align_start, ref1_start, len(ref1))
                print read.query_name, ref_start, cigar, score

        for read in read_bam.fetch(ref, sv_end-mean/2, sv_end+mean/2):
            #print read.query_name, is_poorly_mapped(read)
            if is_poorly_mapped(read):
                align_start, cigar, score =  smapper_core.sw_split(ref1, ref2, read.query_sequence, ref_overlap_len)
                ref_start = local2ref(align_start, ref1_start, len(ref1))
                print read.query_name, ref_start, cigar, score

        '''
        # upstream of DEL
        for read in read_bam.fetch(ref, sv_start-mean, sv_start):
            if read.is_reverse:
                continue
            #print(read.query_name, read.reference_name, read.reference_start, ref_overlap_len)
            mate_seq = read.get_tag("R2")
            mate_seq = utils.rc_seq(mate_seq)
            sv_start, cigar, score =  smapper_core.sw_split(ref1, ref2, mate_seq, ref_overlap_len)
            print read.query_name, cigar, score, "+"
        # downstream of DEL
        for read in read_bam.fetch(ref, sv_end, sv_end+mean):
            if not read.is_reverse:
                continue
            #print(read.query_name, read.reference_name, read.reference_start, ref_overlap_len)
            mate_seq = read.get_tag("R2")
            sv_start, cigar, score =  smapper_core.sw_split(ref1, ref2, mate_seq, ref_overlap_len)
            print read.query_name, cigar, score, "-"
        '''
