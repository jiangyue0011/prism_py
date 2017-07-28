#!/usr/bin/env python

import sys
import pysam
import prism_hash2 as prism
import common

class Cluster(object):
	def __init__(self, line):
		tmp = line.split()
		self.chr1 = tmp[0]
		self.pos1 = int(tmp[1])
		self.chr2 = tmp[2]
		self.pos2 = int(tmp[3])
		self.sv_type = tmp[4]


class SplitReference(object):
	ext_len_up = 100
	ext_len_down = 500
	matrix_len = ext_len_up + ext_len_down
	def __init__(self, cluster, refs):
		self.sv_type = cluster.sv_type
		self.chr1 = cluster.chr1
		self.chr2 = cluster.chr2
		if cluster.sv_type == "DEL":
			self.start1 = cluster.pos1 - self.ext_len_up
			self.end1 = cluster.pos1 + self.ext_len_down
			self.start2 = cluster.pos2 - self.ext_len_down
			self.end2 = cluster.pos2 + self.ext_len_up
			self.overlap_len = self.start2 - self.end1
		elif cluster.sv_type == "DUP":
			self.start1 = cluster.pos2 - self.ext_len_up
			self.end1 = cluster.pos2 + self.ext_len_down
			self.start2 = cluster.pos1 - self.ext_len_down
			self.end2 = cluster.pos1 + self.ext_len_up
			self.overlap_len = self.start1 - self.end2
		else:
			#TODO other types of SVs
			pass
		self.seq1 = refs.fetch(self.chr1, self.start1, self.end1).upper()
		self.seq2 = refs.fetch(self.chr2, self.start2, self.end2).upper()
		self.kmer_dict = common.build_hashtable("N"+self.seq1+self.seq2, 17)
		self.kmer_dict = []


bam_file = sys.argv[1]
ref_file = sys.argv[2]
cluster_file = sys.argv[3]

refs = pysam.FastaFile(ref_file)
bam = pysam.AlignmentFile(bam_file, 'rb')
clusters = []
for line in open(cluster_file, 'r'):
	cluster = Cluster(line)
	if cluster.sv_type != "DUP":
		continue
	print "cluster:", cluster.chr1, cluster.pos1, cluster.chr2, cluster.pos2, cluster.sv_type
	split_ref = SplitReference(cluster, refs)
	print "start,end:", split_ref.start1, split_ref.end1, split_ref.start2, split_ref.end2
	print split_ref.seq1
	print split_ref.seq2

	#split_ref.seq1 = refs.fetch(split_ref.chr1, split_ref.start1, split_ref.end1).upper()
	#split_ref.seq2 = refs.fetch(split_ref.chr2, split_ref.start2, split_ref.end2).upper()
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
		#if read.query_name != "donor_5578_6037_2:0:0_0:0:0_1e38":
		#	continue
		n += 1
		if n % 100 == 0:
			print >> sys.stderr, n
			#break
		print "-------------------------------------------------------"
		print read.tostring(bam)
		split_aln = prism.sw_split(split_ref, read.query_sequence)
		read.cigarstring = split_aln.cigar_string
		read.reference_start = split_aln.mapping_pos-1 # pysam use 0-based coordinate
		read.set_tags([])
		print read.tostring(bam)
		#common.print_sam(split_aln, read)
	sys.exit(1)
	
