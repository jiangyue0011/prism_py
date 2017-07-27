#!/usr/bin/env python
import sys

def message(msg, exit=0):
	print >> sys.stderr, msg
	if exit:
		sys.exit(1)
