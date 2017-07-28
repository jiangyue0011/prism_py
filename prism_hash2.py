#!/usr/bin/env python

import sys
from utils import *
import common
#from numpy import *

FROM_NORTH_NORTH = 0x1
FROM_NORTH_NORTHWEST = 0x2
FROM_WEST_NORTHWEST = 0x3
FROM_WEST_WEST = 0x4
FROM_NORTHWEST_NORTH = 0x5
FROM_NORTHWEST_NORTHWEST = 0x6
FROM_NORTHWEST_WEST = 0x7
FROM_JUMP_NORTHWEST =  0x09
#FROM_NORTHWEST_MAX_NORTH = 0x08
#FROM_NORTHWEST_MAX_WEST = 0xA
FROM_NORTH_JUMP = 0x0B
FROM_WEST_JUMP = 0x0C
FROM_NORTHWEST_JUMP = 0x0D
WORST_SCORE = -10000

class SWCell(object):
	def __init__(self):
		self.score_north = 0
		self.score_west = 0
		self.score_northwest = 0
		self.source_north = 0
		self.source_northwest = 0
		self.source_west = 0
		self.score_jump = WORST_SCORE
		self.source_jump = ""
		self.jump_row_idx = 0
		self.jump_col_idx = 0
		self.score_max = 0
		
		# soucre of next cell
		#self.source_row = 0 
		#self.source_col = 0
		# source of next cell but one
		#self.source_next_row = 0 
		#self.source_next_col = 0


class SplitAlignment(object):
	def __init__(self):
		self.mapping_pos = 0
		self.split_start = 0
		self.split_end = 0
		self.cigars = []

	def build_cigar_string(self, split_ref):
		for idx, cigar in enumerate(self.cigars):
			if cigar.endswith("N"):
				code_len = int(cigar[0:-1])
				print "code", cigar, code_len
				if split_ref.sv_type == "DEL":
					self.cigars[idx] = str(code_len + split_ref.overlap_len) + "N"
				elif split_ref.sv_type == "DUP":
					self.cigars[idx] = str(2*split_ref.matrix_len - code_len + split_ref.overlap_len) + "N"
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


# update each cell's max score and where its source
def update_max_score(max_score, max_from ,cur_score, cur_from):
	if cur_score > max_score:
		return cur_score, cur_from
	else:
		return max_score, max_from

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


def cal_sw_score(sw_matrix, matrix1_max_scores, row_idx, col_idx, match, jump_end):
	# define scores, TODO move to a public clss
	score_gap_open = -47
	score_gap_ext = -7
	if match:
		score_match_mismatch = 10
	else:
		score_match_mismatch = -15
	score_jump = -1
	# calculate scores for each cell
	cell = sw_matrix[row_idx][col_idx]
	cell_north = sw_matrix[row_idx-1][col_idx]
	cell_northwest = sw_matrix[row_idx-1][col_idx-1]
	cell_west = sw_matrix[row_idx][col_idx-1]

	max_score = WORST_SCORE
	max_from = -1
	# calculate score_north, either open a new gap (insertion) or extend previous gap
	score_north_northwest = cell_north.score_northwest + score_gap_open
	max_score, max_from = update_max_score(max_score, max_from, score_north_northwest, FROM_NORTH_NORTHWEST)

	score_north_north = cell_north.score_north + score_gap_ext
	max_score, max_from = update_max_score(max_score, max_from, score_north_north, FROM_NORTH_NORTH)

	score_north_jump = cell_north.score_jump + score_gap_open
	max_score, max_from = update_max_score(max_score, max_from, score_north_north, FROM_NORTH_JUMP)

	cell.score_north = max(score_north_northwest, score_north_north, score_north_jump)
	cell.source_north = max_from
	
	# calculate score_northwest
	score_northwest_north = cell_northwest.score_north + score_match_mismatch
	max_score, max_from = update_max_score(max_score, max_from, score_northwest_north, FROM_NORTHWEST_NORTH)

	score_northwest_northwest = cell_northwest.score_northwest + score_match_mismatch
	max_score, max_from = update_max_score(max_score, max_from, score_northwest_northwest, FROM_NORTHWEST_NORTHWEST)

	score_northwest_west = cell_northwest.score_west + score_match_mismatch
	max_score, max_from = update_max_score(max_score, max_from, score_northwest_west, FROM_NORTHWEST_WEST)

	score_northwest_jump = cell_northwest.score_jump + score_match_mismatch
	max_score, max_from = update_max_score(max_score, max_from, score_northwest_jump, FROM_NORTHWEST_JUMP)

	cell.score_northwest = max(score_northwest_north, score_northwest_northwest, score_northwest_west, score_northwest_jump)
	cell.source_northwest = max_from

	# calculate score_west
	score_west_northwest = cell_west.score_northwest + score_gap_open
	max_score, max_from = update_max_score(max_score, max_from, score_west_northwest, FROM_WEST_NORTHWEST)

	score_west_west = cell_west.score_west + score_gap_ext
	max_score, max_from = update_max_score(max_score, max_from, score_west_west, FROM_WEST_WEST)

	score_west_jump = cell_west.score_northwest + score_gap_open
	max_score, max_from = update_max_score(max_score, max_from, score_west_northwest, FROM_WEST_JUMP)

	cell.score_west = max(score_west_northwest, score_west_west, score_west_jump)
	cell.source_west = max_from

	# calculate score_jump
	if row_idx > 0 and jump_end > 0:
		cell.jump_row_idx = row_idx-1
		cell.jump_col_idx = matrix1_max_scores[row_idx-1][0]
		cell.score_jump = matrix1_max_scores[row_idx-1][1] + score_match_mismatch + score_jump
		# can only jump to a cell which has match/mismatch, makes no sense to jump to a gap
		max_score, max_from = update_max_score(max_score, max_from, cell.score_jump, FROM_JUMP_NORTHWEST)
		cell.source_jump = max_from
		
	cell.score_max = max_score

	# update matrix1_max_scores
	if not jump_end and cell.score_northwest >= matrix1_max_scores[row_idx][1]:
		matrix1_max_scores[row_idx] = [col_idx, cell.score_northwest]

	#print row_idx, col_idx, match, jump_end, cell.source_row, cell.source_col
	#print cell.score_north, cell.score_northwest, cell.score_west, cell.score_jump, cell.score_max

def fill_sw_matrix(sw_matrix, matrix1_max_scores, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="RECTANGLE"):
	diag_len = col_start - row_start		
	for i in range(row_start, row_end):
		if fill_type == "DIAGONAL":
			col_start = i + diag_len
			col_end = col_start + 1
		for j in range(col_start, col_end):
			#if fill_type == "DIAGONAL" and i-row_start != j-col_start:
			#	continue
			match = True if query[i] == db[j] else False
			# jump is only allowed in matrix2
			if j >= matrix1_width:
				jump_end = matrix1_width
			else:
				jump_end = 0
			cal_sw_score(sw_matrix, matrix1_max_scores, i, j, match, jump_end)


# convert split reference alignment to full reference alignment
def build_alignment(split_ref, split_aln):
	chr = split_ref.chr1
	# split_aln.mapping_pos is from pysam.reference_start, this is 0-based
	if split_aln.mapping_pos <= split_ref.matrix_len:
		split_aln.mapping_pos = split_ref.start1 + split_aln.mapping_pos 
	else:
		split_aln.mapping_pos = split_ref.start2 + split_aln.mapping_pos - split_ref.matrix_len
	
	split_aln.cigar_string = split_aln.build_cigar_string(split_ref)
	print "split_aln:", chr, split_aln.mapping_pos, split_aln.cigar_string


# a wrapper which takes SplitReference instance as input
def sw_split(split_ref, read_seq):
	split_aln = sw_split_core(split_ref.seq1, split_ref.seq2, split_ref.kmer_dict, read_seq)
	build_alignment(split_ref, split_aln)
	return split_aln
	

def sw_split_core(ref_seq1, ref_seq2, kmer_dict, read_seq):

	# add 1 "N"	to represent idx 0
	db = "N" + ref_seq1 + ref_seq2
	query = "N" + read_seq
	
	matrix1_width = len(ref_seq1) + 1
	n_rows = len(query)
	n_cols = len(db)

	print "nrow, ncol, mwidth:", n_rows, n_cols, matrix1_width

	# create score matrix
	sw_matrix = [[SWCell() for col in xrange(n_cols)] for row in xrange(n_rows)]
	#initialize the first column
	for i in xrange(n_rows):
		sw_matrix[i][0].source_north = FROM_NORTH_NORTH
		sw_matrix[i][0].score_northwest = WORST_SCORE
		sw_matrix[i][0].score_west = WORST_SCORE
	# this list saves max score for each row in matrix1, this is used when calculating jump score
	matrix1_max_scores = [[0, WORST_SCORE] for i in xrange(n_rows)]
	
	kmer_list = common.query_hashtable(kmer_dict, query, 17)
	row_start = 1
	col_start = 1
	row_end = 1
	col_end = 1
	for kmer in kmer_list:
		row_end = kmer[0]
		col_end = kmer[1]
		# fill the square
		fill_sw_matrix(sw_matrix, matrix1_max_scores, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="RECTANGLE")
		# fill the diagonal
		row_start = row_end		 
		col_start = col_end		 
		row_end = row_start+kmer[2] 
		col_end = col_start+kmer[2]
		fill_sw_matrix(sw_matrix, matrix1_max_scores, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="DIAGONAL")
	# fill the rest area after last kmer
	row_start = row_end		 
	col_start = col_end
	row_end = n_rows 
	col_end = n_cols
	fill_sw_matrix(sw_matrix, matrix1_max_scores, db, query, row_start, row_end, col_start, col_end, matrix1_width, fill_type="RECTANGLE")

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
	for i in range(n_rows):
		for j in range(n_cols):
			cur_cell = sw_matrix[i][j]
			print ",".join([str(x) for x in [cur_cell.score_north, cur_cell.score_northwest, cur_cell.score_west, cur_cell.score_jump, cur_cell.score_max]]),
		print
	'''

	split_aln = backtrace(sw_matrix)
	return split_aln


def backtrace(sw_matrix):
	# starts from the last row (end of read)
	row_idx = len(sw_matrix)-1
	# find the cell with max score
	last_row_scores = [sw_matrix[row_idx][j].score_max for j in xrange(0, len(sw_matrix[0]))]
	max_score = max(last_row_scores)
	col_idx = last_row_scores.index(max_score)
	# find the source where this cell get the max score
	cell = sw_matrix[row_idx][col_idx]
	if max_score == cell.score_north:
		source = cell.source_north
	elif max_score == cell.score_northwest:
		source = cell.source_northwest
	elif max_score == cell.score_west:
		source = cell.source_west
	else:	
		source = cell.source_jump

	split_aln = SplitAlignment()	
	cur_code = "" # cigar code: M, D, I, N
	last_code = ""
	jump_len = 0
	code_len = 0
	code_list = []
	next_cell = cell
	while row_idx > 0:
		cur_cell = next_cell
		print (row_idx, col_idx), cur_cell.score_max, cur_code, last_code, (cur_cell.score_north, cur_cell.score_northwest, cur_cell.score_west, cur_cell.score_jump), cur_cell.jump_row_idx, cur_cell.jump_col_idx, source 
		# cell north
		if (source == FROM_NORTH_NORTH or 
				source == FROM_NORTH_NORTHWEST or
				source == FROM_NORTH_JUMP):
			cur_code = "I"
			row_idx -= 1
		# cell northwest
		elif (source == FROM_NORTHWEST_NORTH or
				source == FROM_NORTHWEST_NORTHWEST or
				source == FROM_NORTHWEST_WEST or
				source == FROM_NORTHWEST_JUMP):
			cur_code = "M"
			row_idx -= 1
			col_idx -= 1
		# cell west
		elif (source == FROM_WEST_NORTHWEST or
				source == FROM_WEST_WEST or
				source == FROM_WEST_JUMP):
			cur_code = "D"
			col_idx -= 1
		# cell jump
		elif (source == FROM_JUMP_NORTHWEST):
			cur_code = "M"
			jump_len = col_idx - cur_cell.jump_col_idx - 1
			split_aln.split_start = cur_cell.jump_col_idx
			split_aln.split_end = col_idx
			row_idx = cur_cell.jump_row_idx
			col_idx = cur_cell.jump_col_idx
		else:
			message("source error:"+str(source), exit=1)
		# generate cigar code
		if cur_code != last_code and last_code != "":
			code_list.append(str(code_len)+last_code)
			code_len = 1
		else:
			code_len += 1
		if source == FROM_JUMP_NORTHWEST:
			code_list.append(str(code_len)+cur_code)
			# add "N"s to cigar
			cur_code = "N"
			code_len = jump_len
		last_code = cur_code

		# continue backtrace
		print "source1:", source
		next_cell = sw_matrix[row_idx][col_idx]
		if (source == FROM_NORTH_NORTH or
				source == FROM_NORTHWEST_NORTH):
			source = next_cell.source_north
		elif (source == FROM_NORTH_NORTHWEST or
				source == FROM_NORTHWEST_NORTHWEST or
				source == FROM_WEST_NORTHWEST or
				source == FROM_JUMP_NORTHWEST):
			source = next_cell.source_northwest
		elif (source == FROM_NORTHWEST_WEST or
					source == FROM_WEST_WEST):
			source = next_cell.source_west
		elif (source == FROM_NORTH_JUMP or
				source == FROM_NORTHWEST_JUMP or
				source == FROM_WEST_JUMP):
			source = next_cell.source_jump
		else:
			message("source error:"+str(source), exit=1)
		print "source2:", source
		#print next_cell.source_north, next_cell.source_northwest, next_cell.source_west, next_cell.source_jump
		
		'''
		if row_idx == cell.source_row + 1 and col_idx == cell.source_col + 1: # match/mismatch
			cur_code = "M"
		elif row_idx == cell.source_row and col_idx == cell.source_col + 1: # deletion
			cur_code = "D"
		elif row_idx == cell.source_row + 1 and col_idx == cell.source_col: # insertion
			cur_code = "I"
		else: # jump
			cur_code = "M"
			is_jump = True
			jump_pos = col_idx
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
			split_aln.split_start = cell.source_col
			split_aln.split_end = col_idx
		last_code = cur_code
		
		print (row_idx, col_idx), cell.score_max, (cell.score_north, cell.score_northwest, cell.score_west), is_jump, cell.source_row, cell.source_col
		#print code_list

		row_idx = cell.source_row
		col_idx = cell.source_col
		'''
	#if cur_code != "M":
	#	message("ERROR: alignment starts with indels.", exit=0)
	code_list.append(str(code_len)+cur_code)
	split_aln.cigars = code_list[::-1]
	split_aln.mapping_pos = col_idx + 1
	split_aln.cigar_string = "".join(code_list[::-1]) 
	print "POS,CIGAR,JUMP:", split_aln.mapping_pos, split_aln.cigar_string, split_aln.cigars, split_aln.split_start, split_aln.split_end
	return split_aln

if __name__ == '__main__':
	#ref_file = sys.argv[1]
	#read_file = sys.argv[2]
	# load all reference sequences into memory. Use indexed refrerence with pysam methods is slower.
	#ref = Reference(ref_file)

	#sw_split(ref_seq1, ref_seq2, read_seq)
	ref1 = "BACCABBBBB"
	ref2 = "BBBDACDFFCAB"
	read = "DACDACCA"
	#ref1 = "BACBBBBB"
	#ref2 = "DDDDAEB"
	#read = "ACAE"
	print ref1, ref2
	print read
	#hash("N"+ref1+ref2, "N"+read, 3)
	split_aln = sw_split_core(ref1, ref2, [], read)

	#sw_split("BACBBBBB","BBBDDFFB","ACDDEEFF")
