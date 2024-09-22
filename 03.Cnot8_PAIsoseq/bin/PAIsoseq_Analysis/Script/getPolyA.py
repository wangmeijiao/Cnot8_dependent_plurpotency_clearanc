#!/usr/bin/env python
import os
import sys
import pysam
import random
import regex as re
from itertools import islice

def scoreContinues(seq):
	# Un-weighted scoring
	score = 0
	for i in range( len(seq) - 1):
		if seq[i] != seq[i+1]:
			score += 1
	return score

if len(sys.argv) != 3:
	msg = ['python',sys.argv[0],'input','window_size']
	msg = ' '.join([str(x) for x in msg])
	print(msg)
	print('\tinput: Input badPolyA file')
	print('\twindow_size: Window size')
	print('\toutput: Stand out')
	quit()

infile = sys.argv[1]
window = int(sys.argv[2])

with open(infile,'r') as fin:
	for line in fin:
		(ID,three_soft_clip,three_soft_clip_35,ccs_full) = line.strip().split("\t")

		# Draft polyA tail from 5'-side of three_soft_clip
		polyA = three_soft_clip 
		window= 10
		window_score = ""
		polyA = polyA + 'A'*10
		# Scan from 5'->3' with slide windows and mark each window. 
		for idx in range( len(polyA) - window + 1):
			window_seq = polyA[idx:(idx + window)]
			window_A = window_seq.count('A')
			window_T = window_seq.count('T')
			window_C = window_seq.count('C')
			window_G = window_seq.count('G')

			score = scoreContinues(window_seq)
			if score <= 4 and window_A + window_T >= 7: 
				window_score = window_score + "0"
			else:
				window_score = window_score + "1"

		#Get the maximum counts of continuous "1" and 
		# the position of last continuous "1" (>=5)
		window_pattern = "1" * int(window/2)
		start_idx = window_score.rfind(window_pattern)
		#If there's just 1 continuous "1"(>=5), and located in the last 2/3 of the polyA
		# (without extending),do not filter out the bad window
		if start_idx == -1:
			polyA = polyA
		else:
			if len(re.findall(window_pattern, window_score)) == 1 and \
			start_idx in range(int((1/3)*(len(window_score)-1)),len(window_score)-1):
				polyA = polyA
		#Else: extract polyA from the next position of the end of last continuous "1"(>=5)
			else:
				polyA = polyA[ (start_idx + len(window_pattern)) :]

		# remove the pseudo 10A
		polyA = polyA[:-10]

		# Filter out the window from the 5' end of polyA until 
		# there's no G or C in first 5nt of polyA
		while( 'C' in polyA[:5] or 'G' in polyA[:5]):
			polyA = polyA[1:]

		# if no polyA, drop the read
		if polyA == '':
			sys.stderr.write(ID + '\n')
			continue
		
		seq1 = three_soft_clip_35[-len(polyA)-35:-len(polyA)] + polyA[:10]
		seq2 = three_soft_clip_35[-len(polyA)-35:-len(polyA)] + polyA[:30]
		seq3 = three_soft_clip_35[-len(polyA)-35:]
		out = [ID, seq1,seq2,seq3, polyA]
		out = '\t'.join([str(x) for x in out])
		print(out)
fin.close()
