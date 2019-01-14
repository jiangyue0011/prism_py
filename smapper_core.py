#!/usr/bin/env python

import sys
import utils
import numpy as np

class SWCell(object):
    def __init__(self):
        self.score_up = 0
        self.score_left = 0
        self.score_diag = 0
        self.score_jump = 0
        self.score_max = 0
        
        self.source_row = 0 
        self.source_col = 0
        # initialize firrst column
        #if col_idx == 0:
        #    self.source_row = row_idx-1 


class Reference(object):
    """
    A dictionary to save whole reference sequences in the ram.
    """
    def __init__(self, ref_file):
        self.ref_dict = {}  # key = chr name; value = sequence
        self.load_reference(ref_file)

    def load_reference(self, ref_file):
        seqs = []
        with open(ref_file, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if seqs:
                        self.ref_dict.setdefault(chrom, "".join(seqs))  
                        seqs = []
                    chrom = line.split()[0][1:]
                else:
                    seqs.append(line.strip())
        self.ref_dict.setdefault(chrom, "".join(seqs))  # last chr 

        utils.message("Reference loaded.")


def hash_kmers(db, query, kmer_len, matrix1_width):
    """
    Hash db kmers to quickly find matched kmers between query and db.
    """
    kmer_dict = {}
    matched_kmers = []
    for i in xrange(1, matrix1_width-kmer_len+1):  # hash kmer in ref1
        kmer_dict.setdefault(db[i:i+kmer_len], []).append(i)
        #print "ref1 kmer", i, db[i:i+kmer_len]
    for i in xrange(matrix1_width, len(db)-kmer_len+1):  # hash kmer in ref2
        kmer_dict.setdefault(db[i:i+kmer_len], []).append(i)
        #print "ref2 kmer", i, db[i:i+kmer_len]
    for i in xrange(1, len(query)-kmer_len+1):  # look for query kmers in db kmer dict
        q_kmer = query[i:i+kmer_len]
        #print "read kmer", i, q_kmer
        if q_kmer in kmer_dict:  # look for query kmer in hasded db
            if len(kmer_dict[q_kmer]) == 1:  # only use unique kmers
                j = kmer_dict[q_kmer][0]
                matched_kmers.append([i, j, 0])  # [query_pos, db_pos, used_for_merge]
    #print matched_kmers
    merged_kmers = []  # save all kmers shared by query and db
    for i, kmer in enumerate(matched_kmers):
        if kmer[2]:
            continue 
        matched_kmers[i][2] = 1
        segment_start_row = matched_kmers[i][0]
        segment_start_col = matched_kmers[i][1]
        segment_len = 1
        #print "    matched i,", i, matched_kmers[i]
        for j in xrange(i+1, len(matched_kmers)):
            # matched_kmers was built in a order of query pos
            if matched_kmers[j][0] - matched_kmers[i][0] == matched_kmers[j][1] - matched_kmers[i][1] == segment_len:
                #print "    matched j,", j, matched_kmers[j]
                segment_len += 1
                matched_kmers[j][2] = 1
            elif matched_kmers[j][0] - matched_kmers[i][0] > segment_len: 
                break
        # dump merged kmer
        merged_kmers.append([segment_start_row, segment_start_col, segment_len+kmer_len-1])  #[row_start, col_start, merged kmer length]
    selected_kmers = []
    for kmer in sorted(merged_kmers, key=lambda x:x[2], reverse=True):
        overlap = False
        for selected_kmer in selected_kmers:
            #if max(selected_kmer[0], kmer[0]) < min(selected_kmer[0]+selected_kmer[2], kmer[0]+kmer[2]) or max(selected_kmer[1], kmer[1]) < min(selected_kmer[1]+selected_kmer[2], kmer[1]+kmer[2]):
            if not ((kmer[0]+kmer[2]-1 < selected_kmer[0] and kmer[1]+kmer[2]-1 < selected_kmer[1]) or 
                    (kmer[0] > selected_kmer[0]+selected_kmer[2]-1 and kmer[1] > selected_kmer[1]+selected_kmer[2]-1)):  # new kmer must locate in the lefttop or rightbottom area of selected kmers
                overlap = True
                break
        if not overlap:
            selected_kmers.append(kmer)
    #print "selected:", selected_kmers
    return sorted(selected_kmers, key=lambda x:x[0])


def cal_sw_score(sw_matrix, row_idx, col_idx, match, jump_end, matrix1_width):
    score_gap_open = -40
    score_gap_ext = -7
    if match:
        score_match_mismatch = 10
    else:
        score_match_mismatch = -15
    score_jump = -1

    cell = sw_matrix[row_idx][col_idx]
    if cell['filled'] == 1:  # the cell's scores have been calculated
        return

    if col_idx != matrix1_width:  # the first column in matrix2 can only have jump score
        cell_up = sw_matrix[row_idx-1][col_idx]
        cell_diag = sw_matrix[row_idx-1][col_idx-1]
        cell_left = sw_matrix[row_idx][col_idx-1]
        cell['score_up_diag'] = cell_up['score_diag'] + score_gap_open + score_gap_ext
        cell['score_up_up'] = max(cell_up['score_up_up'], cell_up['score_up_diag']) + score_gap_ext
        cell['score_diag'] = cell_diag['score_max'] + score_match_mismatch
        cell['score_left_diag'] = cell_left['score_diag'] + score_gap_open + score_gap_ext
        cell['score_left_left'] = max(cell_left['score_left_left'], cell_left['score_left_diag']) + score_gap_ext
    
    jump_idx = 0
    if row_idx > 1 and jump_end > 0:  # calculate jump score
        score_list = [sw_matrix[row_idx-1][j]['score_diag'] for j in range(1, jump_end+1)]
        cell['score_jump'] = max(score_list)
        jump_idx = score_list.index(cell['score_jump'])+1
        #print "jump end:", (row_idx, col_idx), jump_end, jump_idx, cell['score_jump']
        if sw_matrix[row_idx-1][jump_idx]['filled'] == 1:
            cell['score_jump'] += score_jump + score_match_mismatch
            cell['jump_col'] = jump_idx
            #cell['score_jump'] = -999
    #print "cell", row_idx, col_idx, cell
    score_list = [cell['score_diag'], cell['score_up_up'], cell['score_up_diag'], cell['score_left_left'], cell['score_left_diag'], cell['score_jump']]  # "jump" is the last choice
    #print "cell info", row_idx, col_idx, score_list
    cell['score_max'] = max(score_list)
    cell['filled'] = 1

def fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, ref_overlap, fill_type="RECTANGLE"):
    for i in range(row_start, row_end+1):
        for j in range(col_start, col_end+1):
            if fill_type == "DIAGONAL" and i-row_start != j-col_start:
                continue
            match = True if query[i] == db[j] else False
            if j >= matrix1_width:  # jump is only allowed in matrix2
                jump_end = min(j - ref_overlap -1, matrix1_width)
            else:
                jump_end = 0
            score = cal_sw_score(sw_matrix, i, j, match, jump_end, matrix1_width)
    

def sw_split(ref_seq1, ref_seq2, read_seq, ref_overlap):

    # add 1 "N"  to represent idx 0
    db = "N" + ref_seq1 + ref_seq2
    query = "N" + read_seq
    matrix1_width = len(ref_seq1) + 1
    n_rows = len(query)
    n_cols = len(db)
    swcell = np.dtype({'names':['score_up_diag', 'score_up_up', 'score_left_diag', 'score_left_left', 'score_diag', 'score_jump', 'score_max', 'jump_col', 'filled'], 
                       'formats':['i', 'i', 'i', 'i', 'i', 'i', 'i', 'i', 'i']}, align = True)
    sw_matrix = np.full((n_rows, n_cols), -999, dtype=swcell)

    '''
    # Test the speed of ndarray and tuple list
    n_rows = 10000
    n_cols = 10000
    if 1:
        swcell = np.dtype({'names':['score_up_diag', 'score_up_up', 'score_left_diag', 'score_left_left', 'score_diag', 'score_jump', 'score_max', 'jump_col', 'filled'], 
                           'formats':['i', 'i', 'i', 'i', 'i', 'i', 'i', 'i', 'i']}, align = True)
        sw_matrix = np.full((n_rows, n_cols), -999, dtype=swcell)
        for i in range(10):
            print sw_matrix[i][0]['score_diag']
    else:
        sw_matrix = [[(1,2,3,4,5)]*n_rows for i in range(n_cols)]
        for i in range(10):
            print sw_matrix[i][0][2]


    sys.exit(1)
    '''

    # initialize first row and first column
    for i in xrange(0, n_rows):
        sw_matrix[i][0]['score_max'] = -(40 + i*7)
    for j in xrange(0, n_cols):
        sw_matrix[0][j]['score_max'] = 0

    #kmer_list = []
    kmer_list = hash_kmers(db, query, 19, matrix1_width)
    #print kmer_list
    row_start = 1
    col_start = 1
    #print_sw_matrix(sw_matrix)
    for kmer in kmer_list:
        row_end = kmer[0]
        col_end = kmer[1]
        #print "kmer1", kmer, row_start, row_end, col_start, col_end
        # fill the square
        fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, ref_overlap, fill_type="RECTANGLE")
        # fill the diagonal
        row_start = kmer[0]       
        col_start = kmer[1]      
        row_end = row_start + kmer[2] - 1
        col_end = col_start + kmer[2] - 1
        #print "kmer2", row_start, row_end, col_start, col_end
        fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, ref_overlap, fill_type="DIAGONAL")
        row_start = row_end        
        col_start = col_end        
    # fill the rest area after last kmer
    row_end = n_rows - 1
    col_end = n_cols - 1
    #print "kmer3", row_start, row_end, col_start, col_end
    fill_sw_matrix(sw_matrix, db, query, row_start, row_end, col_start, col_end, matrix1_width, ref_overlap, fill_type="RECTANGLE")
    #print_sw_matrix(sw_matrix)

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
    align_start, cigar, score = backtrace(sw_matrix)
    print_alignment(db, query, align_start, cigar)

    return align_start, cigar, score


def backtrace(sw_matrix):
    """
    Backtrace of the sw-matirx. It starts from the first max score cell of the last row. When the max score
    if from diag cell, the score can be any of the diag cell's scores. When the max score is score_up_diag,
    it means this score is from up cell, and the score of up cell is from up cell's diag cell. In this case,
    the next two steps are fixed. The path has to go up then go diag.  
    """
    row_idx = len(sw_matrix)-1
    last_row_scores = [sw_matrix[row_idx][j]['score_max'] for j in xrange(0, len(sw_matrix[0]))]
    max_score = max(last_row_scores)
    col_idx = last_row_scores.index(max_score)
    #print "max score cell:", row_idx, col_idx, max_score
    
    cur_code = ""  # cigar code: M, D, I, N
    last_code = ""
    code_len = 0
    code_list = []
    next_cell = ""  # to force next cell go according to left_left left_diag up_up up_diag
    jump_len = 0

    while row_idx > 0:
        if row_idx == 1:
            align_pos = col_idx
        #print "backtrace", (row_idx, col_idx), next_cell
        cell = sw_matrix[row_idx][col_idx]
        if not next_cell:
            if cell['score_max'] == cell['score_diag']:  # match/mismatch
                row_idx -= 1
                col_idx -= 1
                cur_code = "M"
            elif cell['score_max'] == cell['score_left_left']:  # deletion
                col_idx -= 1
                cur_code = "D"
                next_cell = "left"
            elif cell['score_max'] == cell['score_left_diag']:  # deletion
                col_idx -= 1
                cur_code = "D"
                next_cell = "diag"
            elif cell['score_max'] == cell['score_up_up']:  # insertion
                row_idx -= 1
                cur_code = "I"
                next_cell = "up"
            elif cell['score_max'] == cell['score_up_diag']:  # insertion
                row_idx -= 1
                cur_code = "I"
                next_cell = "diag"
            else:  # jump
                cur_code = "M"
                jump_len = col_idx - cell['jump_col'] - 1
                row_idx -= 1
                col_idx = cell['jump_col']
                next_cell = "diag"
        elif next_cell == "diag":
            row_idx -= 1
            col_idx -= 1
            cur_code = "M"
            next_cell = ""
        elif next_cell == "left":
            col_idx -= 1
            cur_code = "D"
            if cell['score_max'] == cell['score_left_diag']:
                next_cell = "diag"
        elif next_cell == "up":
            row_idx -= 1
            cur_code = "I"
            if cell['score_max'] == cell['score_up_diag']:
                next_cell = "diag"
    
        # build CIGAR string
        if cur_code != last_code and last_code != "":
            code_list.append(str(code_len)+last_code)
            code_len = 1
        else:
            code_len += 1
        if jump_len: # for the case that ref1 and ref2 are continuous
            code_list.append(str(code_len)+cur_code)
            cur_code = "N"
            code_len = jump_len
            jump_len = 0
        last_code = cur_code
    code_list.append(str(code_len)+cur_code)
    cigar = "".join(code_list[::-1])
    if cur_code != "M":
        utils.message("ERROR: alignment starts with indels." + cur_code, exit=1)

    return align_pos, cigar, max_score


def print_sw_matrix(sw_matrix):
    """
    Print socre_max for each cell in the format of a matrix
    """
    for i in xrange(0, sw_matrix.shape[0]):
        for j in xrange(0, sw_matrix.shape[1]):
            print str(sw_matrix[i][j]['score_max']).rjust(4),
        print


def print_alignment(db, query, pos, cigar):
    """
    Print alignment of db and query
    """
    #print db, query, pos, cigar
    out_query = ""
    out_db = ""
    db_pos = pos
    query_pos = 1

    cigar_len = ""
    for s in cigar:
        if s.isdigit():
            cigar_len += s
        else:
            cigar_len = int(cigar_len)
            if s == "M":
                out_db += db[db_pos:db_pos+cigar_len]
                db_pos += cigar_len
                out_query += query[query_pos:query_pos+cigar_len]
                query_pos += cigar_len
            elif s == "N" or s == "D":
                if cigar_len <= 10:
                    out_db += db[db_pos:db_pos+cigar_len]
                    db_pos += cigar_len
                    if s == "N":
                        out_query += "=" * cigar_len
                    else:
                        out_query += "_" * cigar_len
                else:
                    out_db += db[db_pos:db_pos+5] + "..." + db[db_pos+cigar_len-5:db_pos+cigar_len]
                    db_pos += cigar_len
                    if s == "N":
                        out_query += "=" * (10 + len("..."))
                    else:
                        out_query += "_" * (10 + len("..."))
            elif s == "I":
                out_db += "_" * cigar_len
                out_query += query[query_pos:query_pos+cigar_len]
                query_pos += cigar_len
                
            else:
                raise Exception("Unkonwn cigar code:", s)
            cigar_len = ""
    print out_db
    print out_query
                

if __name__ == '__main__':
    ref_file = sys.argv[1]
    #read_file = sys.argv[2]
    # load all reference sequences into memory. Use indexed refrerence with pysam methods is slower.
    ref = Reference(ref_file)
    print ref.ref_dict

    #sw_split(ref_seq1, ref_seq2, read_seq)
    ref1 = "BACCABBB" *1 
    ref2 = "BBBDACDFFCAB" * 1
    read = "ACCADACDEEFFCA"
    #ref1 = "BBBBAAAA"
    #ref2 = "BBBBAAAA"
    ref_overlap = 0
    #read = "AAABBB"
    #ref1 = "BACB"
    #ref2 = "BFFB"
    #read = "ACFF"
    print ref1, ref2
    print read
    print len(ref1), len(ref2), len(read)
    for i in xrange(0, 100):
        sw_split(ref1, ref2, read, ref_overlap)
