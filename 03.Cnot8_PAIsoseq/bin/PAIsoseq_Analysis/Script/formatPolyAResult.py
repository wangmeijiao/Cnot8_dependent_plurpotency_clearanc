#!/usr/bin/env python
import os
import sys

if len(sys.argv) != 3:
	msg = ['python',sys.argv[0], 'sample.polyA','sample.bam.featureCounts','>sample.polyA_tail_length.txt']
	msg = ' '.join([ str(x) for x in msg])
	print(msg)
	quit()

infile = sys.argv[1]
fcounts = sys.argv[2]

fc_dict = {}
with open(fcounts, 'r') as fc:
	for line in fc:
		(name, status, strand, gene) = line.strip('\n').split('\t')
		if status == 'Assigned':
			fc_dict[name] = gene
		else:
			fc_dict[name] = 'Unknown'

with open(infile, 'r') as infile:
	for line in infile:
		(name, seq) = line.strip('\n').split('\t')
		A = seq.count('A')
		T = seq.count('T')
		C = seq.count('C')
		G = seq.count('G')
		nonA = T + C + G
		(read_name, sample, read_pass) = name.split(':')

		gene = ''
		if name in fc_dict:
			gene = fc_dict[name]
		else:
			print('Error!')
			quit()

		out = [sample, name, gene, read_pass, 1, A, T, C, G, nonA, 0, seq]
		out = '\t'.join([str(x) for x in out])
		print(out)
