#!/usr/bin/env python

import sys
class SWCell(object):
	def __init__(self):
		self.score_up = 0
		self.score_left = 0
		self.score_diag = 0
		self.score_jump = -10000
		self.score_max = 0
		
		self.source_row = 0 
		self.source_col = 0
		# initialize firrst column
		#if col_idx == 0:
		#	self.source_row = row_idx-1 


class SplitAlignment(object):
	def __init__(self):
		self.mapping_pos = 0
		self.split_start = 0
		self.split_end = 0
		self.cigars = []
		self.cigar_string = ""

	def build_cigar_string(self, overlap_len):
		for idx, cigar in enumerate(self.cigars):
			if cigar.endswith("N"):
				code_len = int(cigar[0:-1])
				print "code", cigar, code_len
				self.cigars[idx] = str(code_len + overlap_len) + "N"
				
				break
		return "".join(self.cigars)


class Reference(object):
	def __init__(self, ref_file):
		self.chr_list = []
		self.chr_bp_list = []
		self.ref_seq = ""
		self.load_reference(ref_file)
	# ref_file is a fasta
	def load_reference(self, ref_file):
		refs = []
		ref_len = 0
		chr_len = 0
		with open(ref_file, "r") as f:
			for line in f:
				if line.startswith(">"):
					if chr_len:
						#self.chr_map.setdefault(chr_name, [ref_len, ref_len+chr_len])
						self.chr_list.append(chr_name)
						self.chr_bp_list.append(ref_len+chr_len)
					chr_name = line.split()[0][1:]
					ref_len += chr_len
					chr_len = 0
				else:
					refs.append(line.strip())
					chr_len += len(line.strip())
		if chr_len:
			self.chr_list.append(chr_name)
			self.chr_bp_list.append(ref_len+chr_len)

		self.ref_seq = "".join(refs)
		message("Reference loaded.")
		print self.chr_list, self.chr_bp_list
		print self.ref_seq

# hash split reference
def build_hashtable(db, kmer_len):
	kmer_dict = {}
	for i in xrange(0, len(db)-kmer_len):
		kmer_dict.setdefault(db[i:i+kmer_len], []).append(i)
	return kmer_dict


def query_hashtable(kmer_dict, query, kmer_len):
	# query kmer which are found in db kmer_dict
	matched_kmers = []
	for i in xrange(0, len(query)-kmer_len):
		kmer = query[i:i+kmer_len]
		#print i, kmer
		if kmer in kmer_dict:
			for j in sorted(kmer_dict[kmer]):
				matched_kmers.append([i, j, 0])
		
	merged_kmers = []
	for i, kmer in enumerate(matched_kmers):
		if kmer[2]:
			continue 
		matched_kmers[i][2] = 1
		segment_start_row =  matched_kmers[i][0]
		segment_start_col =  matched_kmers[i][1]
		segment_len = 1
		for j in xrange(i+1, len(matched_kmers)):
			if matched_kmers[j][0] - matched_kmers[i][0] == matched_kmers[j][1] - matched_kmers[i][1] == segment_len:
				print "    matched j,", j, matched_kmers[j]
				segment_len += 1
				matched_kmers[j][2] = 1
			elif matched_kmers[j][0] - matched_kmers[i][0] > segment_len:
				break
		# dump merged kmer
		merged_kmers.append([segment_start_row, segment_start_col, segment_len+kmer_len-1])
	selected_kmers = []
	for kmer in sorted(merged_kmers, key=lambda x:x[2], reverse=True):
		overlap = False
		for selected_kmer in selected_kmers:
			if max(selected_kmer[0], kmer[0]) < min(selected_kmer[0]+selected_kmer[2], kmer[0]+kmer[2]) or max(selected_kmer[1], kmer[1]) < min(selected_kmer[1]+selected_kmer[2], kmer[1]+kmer[2]):
				print "overlap:", selected_kmer
				overlap = True
				break
		if not overlap:
			selected_kmers.append(kmer)
	return selected_kmers

def print_sam(split_aln, read):
	print "\t".join([str(x) for x in [read.query_name, read.flag, read.reference_name, split_aln.mapping_pos, split_aln.cigar_string, read.next_reference_name, read.next_reference_start, read.query_sequence, read.query_qualities]])
	
