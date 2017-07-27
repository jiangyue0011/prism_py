#!/usr/bin/env python

import sys
import utils

#from numpy import *

class SWCell(object):
	def __init__(self, row_idx, col_idx):
		self.score_up = 0
		self.score_left = 0
		self.score_diag = 0
		self.score_jump = 0
		self.score_max = 0
		
		self.source_row = 0 
		self.source_col = 0
		# initialize firrst column
		if col_idx == 0:
			self.source_row = row_idx-1 


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


def hash(db, query, kmer_len):
	kmer_dict = {}
	matched_kmers = []
	for i in xrange(0, len(db)-kmer_len):
		kmer_dict.setdefault(db[i:i+kmer_len], []).append(i)
	'''
	for i in xrange(0, len(ref1)-kmer_len):
		kmer_dict.setdefault(ref1[i:i+kmer_len], []).append(i)
	for i in xrange(0, len(ref2)-kmer_len):
		kmer_dict.setdefault(ref2[i:i+kmer_len], []).append(i+len(ref1))
	'''
	for i in xrange(0, len(query)-kmer_len):
		kmer = query[i:i+kmer_len]
		print i, kmer
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


def cal_sw_score(sw_matrix, row_idx, col_idx, match, jump_end):
	score_gap_open = -40
	score_gap_ext = -7
	if match:
		score_match_mismatch = 10
	else:
		score_match_mismatch = -15
	score_jump = -1

	cell = sw_matrix[row_idx][col_idx]
	cell_up = sw_matrix[row_idx-1][col_idx]
	cell_diag = sw_matrix[row_idx-1][col_idx-1]
	cell_left = sw_matrix[row_idx][col_idx-1]

	cell.score_up = max(cell_up.score_diag + score_gap_open, cell_up.score_up + score_gap_ext)
	cell.score_diag = cell_diag.score_max + score_match_mismatch
	cell.score_left = max(cell_left.score_diag + score_gap_open, cell_left.score_left + score_gap_ext)

	jump_idx = 0
	if row_idx > 0 and jump_end > 0:
		score_list = [sw_matrix[row_idx-1][j].score_diag for j in range(1, jump_end)]
		cell.score_jump = max(score_list)
		jump_idx = score_list.index(cell.score_jump)+1
		cell.score_jump += score_match_mismatch + score_jump
	score_list = [cell.score_diag, cell.score_jump, cell.score_up, cell.score_left]
	cell.score_max = max(score_list)
	max_idx = score_list.index(cell.score_max)
	if max_idx == 0: # max score from diag
		cell.source_row = row_idx-1
		cell.source_col = col_idx-1
	elif max_idx == 1: # max score from jump
		cell.source_row = row_idx-1 
		cell.source_col = jump_idx
	elif max_idx == 2: # max score from up
		cell.source_row = row_idx-1
		cell.source_col = col_idx
	else: # max score from left
		cell.source_row = row_idx
		cell.source_col = col_idx-1

	#print row_idx, col_idx, match, jump_end, cell.source_row, cell.source_col
	#print cell.score_up, cell.score_diag, cell.score_left, cell.score_jump, cell.score_max


def fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="RECTANGLE"):
	for i in range(row_start, row_end):
		for j in range(col_start, col_end):
			if fill_type == "DIAGONAL" and i-row_start != j-col_start:
				continue
			match = True if query[i] == db[j] else False
			# jump is only allowed in matrix2
			if j >= matrix1_width:
				jump_end = matrix1_width
			else:
				jump_end = 0			
			score = cal_sw_score(sw_matrix, i, j, match, jump_end)


def sw_split(split_ref, read_seq):
	return sw_split_core(split_ref.seq1, split_ref.seq2, read_seq)
	

def sw_split_core(ref_seq1, ref_seq2, read_seq):

	# add 1 "N"	to represent idx 0
	db = "N" + ref_seq1 + ref_seq2
	query = "N" + read_seq
	
	matrix1_width = len(ref_seq1) + 1
	n_rows = len(query)
	n_cols = len(db)

	print "nrow, ncol, mwidth:", n_rows, n_cols, matrix1_width

	sw_matrix = [[SWCell(row, col) for col in range(n_cols)] for row in range(n_rows)]
	#kmer_list = []
	kmer_list = hash(db, query, 17)
	row_start = 1
	col_start = 1
	row_end = 1
	col_end = 1
	for kmer in kmer_list:
		row_end = kmer[0]
		col_end = kmer[1]
		# fill the square
		fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="RECTANGLE")
		# fill the diagonal
		row_start = row_end		 
		col_start = col_end		 
		row_end = row_start+kmer[2] 
		col_end = col_start+kmer[2]
		fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="DIAGONAL")
	# fill the rest area after last kmer
	row_start = row_end		 
	col_start = col_end
	row_end = n_rows 
	col_end = n_cols
	fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="RECTANGLE")

	#sw_matrix = zeros((len(read)+1, len(ref1)+len(ref2)+1))
	'''
	for i in range(1, n_rows):
		for j in range(1, n_cols):
			match = True if query[i] == db[j] else False
			# jump is only allowed in matrix2
			if j >= matrix1_width:
				jump_end = matrix1_width
			else:
				jump_end = 0			
			score = cal_sw_score(sw_matrix, i, j, match, jump_end)
	'''
	#for i in range(n_rows):
	#	for j in range(n_cols):
	#		print sw_matrix[i][j].score_max,
	#	print
	backtrace(sw_matrix)


def backtrace(sw_matrix):
	row_idx = len(sw_matrix)-1
	last_row_scores = [sw_matrix[row_idx][j].score_max for j in xrange(0, len(sw_matrix[0]))]
	max_score = max(last_row_scores)
	col_idx = last_row_scores.index(max_score)
	print "max score cell:", row_idx, col_idx, max_score
	
	cur_code = "" # cigar code: M, D, I, N
	last_code = ""
	code_len = 0
	code_list = []
	while row_idx > 0:
		is_jump = False
		cell = sw_matrix[row_idx][col_idx]
		if row_idx == cell.source_row + 1 and col_idx == cell.source_col + 1: # match/mismatch
			cur_code = "M"
		elif row_idx == cell.source_row and col_idx == cell.source_col + 1: # deletion
			cur_code = "D"
		elif row_idx == cell.source_row + 1 and col_idx == cell.source_col: # insertion
			cur_code = "I"
		else: # jump
			cur_code = "M"
			is_jump = True
		print cur_code, last_code,

		if cur_code != last_code and last_code != "":
			code_list.append(str(code_len)+last_code)
			code_len = 1
		else:
			code_len += 1

		if is_jump:
			code_list.append(str(code_len)+cur_code)
			cur_code = "N"
			code_len = col_idx - cell.source_col - 1
			#code_list.append(str(code_len)+cur_code)
		last_code = cur_code
		
		print (row_idx, col_idx), cell.score_max

		row_idx = cell.source_row
		col_idx = cell.source_col
	if cur_code != "M":
		message("ERROR: alignment starts with indels.", exit=0)
	code_list.append(str(code_len)+cur_code)
	print code_list
	print "CIGAR:","".join(code_list[::-1])

if __name__ == '__main__':
	#ref_file = sys.argv[1]
	#read_file = sys.argv[2]
	# load all reference sequences into memory. Use indexed refrerence with pysam methods is slower.
	#ref = Reference(ref_file)

	#sw_split(ref_seq1, ref_seq2, read_seq)
	ref1 = "BACCABBBBB"
	ref2 = "BBBDACDFFCAB"
	read = "ACCADACDEEFFCA"
	print ref1, ref2
	print read
	#hash("N"+ref1+ref2, "N"+read, 3)
	sw_split_core(ref1, ref2, read)

	#sw_split("BACBBBBB","BBBDDFFB","ACDDEEFF")
