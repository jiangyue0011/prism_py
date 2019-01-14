#!/usr/bin/env python
import sys

def message(msg, exit=0):
    print >> sys.stderr, msg
    if exit:
        sys.exit(1)


def rc_seq(seq0):
    seq = seq0.upper()
    tmpseq = seq.replace("A","a")
    seq = tmpseq
    tmpseq = seq.replace("T","A")
    seq = tmpseq
    tmpseq = seq.replace("a","T")
    seq = tmpseq
    tmpseq = seq.replace("G","g")
    seq = tmpseq
    tmpseq = seq.replace("C","G")
    seq = tmpseq
    tmpseq = seq.replace("g","C")
    seq = tmpseq[::-1]
    return seq

