#!/usr/bin/env python
import os
import sys

if len(sys.argv) != 3:
	print('python infile outdir')
	quit()

infile = sys.argv[1]
outdir = sys.argv[2]

last_gene = ''
out = ''
with open(infile, 'r') as f:
	for line in f:
		(name, seq10, seq30, polyA, gene) = line.strip('\n').split('\t')
		if last_gene == '':
			last_gene = gene
		if gene != last_gene:
			outfile = outdir + '/' + last_gene
			of = open(outfile, 'w')
			of.write(out)
			of.close()
			out = '>' + name + '\n' + seq10 + '\n'
		else:
			line_out = '>' + name + '\n' + seq10 + '\n'
			out = out + line_out
		last_gene = gene

outfile = outdir + '/' + last_gene
of = open(outfile, 'w')
of.write(out)
of.close()
