#!/usr/bin/env python

import sys
import pysam
import prism_hash2 as prism
import core

class Cluster(object):
	def __init__(self, line):
		tmp = line.split()
		self.chr1 = tmp[0]
		self.pos1 = int(tmp[1])
		self.chr2 = tmp[2]
		self.pos2 = int(tmp[3])
		self.type = tmp[4]


class SplitReference(object):
	ext_len_up = 100
	ext_len_down = 500
	def __init__(self, cluster):
		self.chr1 = cluster.chr1
		self.chr2 = cluster.chr2
		if cluster.type == "DEL":
			self.start1 = cluster.pos1 - self.ext_len_up
			self.end1 = cluster.pos1 + self.ext_len_down
			self.start2 = cluster.pos2 - self.ext_len_down
			self.end2 = cluster.pos2 + self.ext_len_up
			self.overlap = self.start2 - self.end1
			self.seq1 = ""
			self.seq2 = ""
		else:
			#TODO other types of SVs
			pass
		self.kmer_dict = core.build_hashtable("N"+self.seq1+self.seq2, 17)
		self.kmer_dict = []


bam_file = sys.argv[1]
ref_file = sys.argv[2]
cluster_file = sys.argv[3]

refs = pysam.FastaFile(ref_file)
bam = pysam.AlignmentFile(bam_file, 'rb')
clusters = []
for line in open(cluster_file, 'r'):
	cluster = Cluster(line)
	if cluster.type != "DEL":
		continue
	print "cluster:", cluster.chr1, cluster.pos1, cluster.chr2, cluster.pos2, cluster.type
	split_ref = SplitReference(cluster)
	print "start,end:", split_ref.start1, split_ref.end1, split_ref.start2, split_ref.end2
	split_ref.seq1 = refs.fetch(split_ref.chr1, split_ref.start1, split_ref.end1).upper()
	split_ref.seq2 = refs.fetch(split_ref.chr2, split_ref.start2, split_ref.end2).upper()
	#print split_ref.seq1
	#print split_ref.seq2
	
	#if cluster.type == "DEL":
	#	ref1 = refs.fetch(cluster.chr1, cluster.pos1-100, cluster.pos1+500)
	#	ref2 = refs.fetch(cluster.chr2, cluster.pos2-500, cluster.pos2+100)
	read_list = []
	for read in bam.fetch(cluster.chr1, cluster.pos1-100, cluster.pos1+500):
		if read.cigarstring != str(read.query_length)+"M":
			read_list.append(read)
	for read in bam.fetch(cluster.chr2, cluster.pos2-500, cluster.pos2+100):
		if read.cigarstring != str(read.query_length)+"M":
			read_list.append(read)
	print "#reads:", len(read_list)
	n = 0
	for read in read_list:
		if read.query_name != "donor_4405_4935_2:0:0_1:0:0_697b":
			continue
		n += 1
		if n % 100 == 0:
			print >> sys.stderr, n
			#break
		print "-------------------------------------------------------"
		print read.tostring(bam)
		split_aln = prism.sw_split(split_ref, read.query_sequence)
	#sys.exit(1)
	
