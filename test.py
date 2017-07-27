#!/usr/bin/env python
import sys
import numpy as np
import time
import random
class Cell(object):
	def __init__(self):
		self.score = 0
		self.best = 0

#matrix = np.array([Cell() for i in range(9)])
n_cols = 3
n_rows = 6
#sw_matrix = np.array([[Cell() for col in range(n_cols)] for row in range(n_rows)])
#sw_matrix = [[Cell() for col in range(n_cols)] for row in range(n_rows)]
#sw_matrix = [[Cell() for col in range(n_cols)] for row in range(n_rows)]
#sw_matrix = np.full((n_cols, n_rows), Cell())
swcell = np.dtype({'names':['score_up', 'score_left', 'score_diag', 'score_jump', 'score_max', 'source_row', 'source_col'], 'formats':['i', 'i', 'i', 'i', 'i', 'i', 'i']}, align = True)
sw_matrix = np.zeros(shape=(n_rows, n_cols),dtype=swcell)
#sw_matrix = np.full((n_cols, n_rows), 1)
#print sw_matrix[0][1].score
#print sw_matrix[1][1]
print sw_matrix[1][1][] 

li = [random.random() for i in range(1000000)]
t1 = time.time()
max(li)
t2 = time.time()
print t2-t1

t1 = time.time()
a = (np.array(li))
t2 = time.time()
print t2-t1

t1 = time.time()
np.max(a)
t2 = time.time()
print t2-t1
