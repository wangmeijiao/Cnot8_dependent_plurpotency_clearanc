#!/usr/bin/env python
import pysam
import sys

if len(sys.argv) != 2:
	msg = ['python', sys.argv[0], 'bamfile', '1 >out', ' 2 >err']
	msg = ' '.join([str(x) for x in msg])
	print(msg)
	quit()

def revComp(s):
	revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
	return ''.join([revCDict[x] for x in s])[::-1]

infile = sys.argv[1]
bamfile = pysam.AlignmentFile(infile, 'rb')
for aln in bamfile:
	if aln.is_unmapped:
		sys.stderr.write(aln.query_name + '\tunmapped\n')
		continue

	if aln.is_secondary:
		sys.stderr.write(aln.query_name + '\tsecondary\n')
		continue

	if aln.is_supplementary:
		sys.stderr.write(aln.query_name + '\tsupplementary\n')
		continue

	soft_clip_len = 0
	soft_clip_seq = ''
	if aln.is_reverse == False: # 5'->cDNA->polyA->3'
		last_cigar = aln.cigartuples[-1][0]
		if last_cigar != 4:
			soft_clip_len = 0
		else:
			soft_clip_len = aln.cigartuples[-1][1]
			soft_clip_seq = aln.query_sequence[-1*soft_clip_len:]

	else: # 3'->polyT->cDNA->5'
		first_cigar = aln.cigartuples[0][0]
		if first_cigar != 4:
			soft_clip_len = 0
		else:
			soft_clip_len = aln.cigartuples[0][1]
			soft_clip_seq = aln.query_sequence[:soft_clip_len]
			soft_clip_seq = revComp(soft_clip_seq)

	if soft_clip_seq != 0:
		print('>' + aln.query_name)
		print(soft_clip_seq)
