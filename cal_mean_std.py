#!/usr/bin/env python
import sys
import numpy
import pysam

bam_file = sys.argv[1]
bam = pysam.AlignmentFile(bam_file, "rb")
isize_list = []

iszie_upper = 1000
iszie_lower = 100
sampleing_cnt = 100000
read_cnt = 0
for read in bam.fetch():
	if read.reference_name == read.next_reference_name and read.reference_start < read.next_reference_start and iszie_lower < read.template_length < iszie_upper:
		isize_list.append(read.template_length)
		read_cnt += 1
		if read_cnt == sampleing_cnt:
			break

print >> sys.stderr, read_cnt, "pairs loaded."
print "Mean:\t", numpy.average(isize_list)
print "Stddev:\t",numpy.std(isize_list)





