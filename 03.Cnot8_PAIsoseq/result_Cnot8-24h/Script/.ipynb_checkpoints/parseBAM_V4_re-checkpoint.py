#!/usr/bin/env python
import re
import sys
import pysam

if len(sys.argv) != 2:
	msg = ['python', sys.argv[0], 'bamfile','1 >out', ' 2 >err']
	msg = ' '.join([str(x) for x in msg])
	print(msg)
	print("\tbamfile: Input bam file")
	print('Version: 0.0.1, 2020-05-12')
	print('Author: Hu Nie, niehu@genetics.ac.cn')
	quit()

def revComp(s):
	revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
	return ''.join([revCDict[x] for x in s])[::-1]

# Read arguments
bamFile = sys.argv[1]

# Read bam file, and extract reads
bamfile = pysam.AlignmentFile(bamFile, 'rb')
for aln in bamfile:
	if aln.is_unmapped:
		continue
	if aln.is_secondary:
		continue
	if aln.is_supplementary:
		continue

	soft_clip_len = 0
	soft_clip_seq = ''
	query_seq = ''
	if aln.is_reverse == False : # 5'->cDNA->polyA->3'
		query_seq = aln.query_sequence
		last_cigar = aln.cigartuples[-1][0]
		if last_cigar != 4:
			soft_clip_len = 0
		else:
			soft_clip_len = aln.cigartuples[-1][1]
			soft_clip_seq = aln.query_sequence[-1*soft_clip_len:]
	elif aln.is_reverse == True : # 5'->polyT->cDNA->3'
		query_seq = revComp(aln.query_sequence)
		first_cigar = aln.cigartuples[0][0]
		if first_cigar != 4:
			soft_clip_len = 0
		else:
			soft_clip_len = aln.cigartuples[0][1]
			soft_clip_seq = aln.query_sequence[:soft_clip_len]
			soft_clip_seq = revComp(soft_clip_seq)

	# copy 3'-soft clip
	old_soft_clip_seq = soft_clip_seq

	# make sure the first 5 bases in 3'-soft clip do not contain C or G
	while 'C' in soft_clip_seq[:5] or 'G' in soft_clip_seq[:5]:
		soft_clip_seq = soft_clip_seq[1:]

	# extend polyA
	for i in range( len(query_seq) - len(soft_clip_seq) )[::-1]:
		if query_seq[i] == 'A':
			soft_clip_seq = 'A' + soft_clip_seq
		else:
			break

	out = []
	thres = 3
	window_size = 10

	# calculate gc by window
	if len(soft_clip_seq) <= 10: # 3'-soft clip length <= 10
		if soft_clip_seq.count('C') + soft_clip_seq.count('G') > thres:
			out = [ "FALSE", aln.query_name, "New:", "Old:" + old_soft_clip_seq]
		else:
			out = [ "TRUE", aln.query_name, "New:" + soft_clip_seq, "Old:" + old_soft_clip_seq]
	else: # 3'-soft clip length > 10
		window_gc_list = []
		window_gc_str = ""
		for i in range( len(soft_clip_seq) - window_size + 1):
			seq = soft_clip_seq[i:i+window_size]
			gc = seq.count('C') + seq.count('G')
			window_gc_list.append(gc)

			if gc > thres:
				window_gc_str = window_gc_str + "F"
			else:
				window_gc_str = window_gc_str + "T"

		# compress window gc str
		zip_window_gc_str = ""
		for i in range( len(window_gc_str) ):
			if zip_window_gc_str != "":
				if window_gc_str[i] == re.split("\d+", zip_window_gc_str)[-2]:
					continue
			count = 0
			for j in range(i, len(window_gc_str)):
				if window_gc_str[i] == window_gc_str[j]:
					count = count + 1
				else:
					break
			zip_window_gc_str = zip_window_gc_str + window_gc_str[i] + str(count)

		# find pattern
		pattern = re.sub("\d+", "" ,zip_window_gc_str)
		if pattern == "F":
			out = ["FALSE", aln.query_name, "New:", "Old:" + old_soft_clip_seq]
		elif pattern == "T":
			out = ["TRUE", aln.query_name, "New:" + soft_clip_seq, "Old:" + old_soft_clip_seq]
		elif pattern == "FT":
			flag = 0
			for i in range( len(window_gc_list) ):
				if window_gc_list[i] == 0:
					flag = 1
					soft_clip_seq = soft_clip_seq[i:]
					out = ["TRUE", aln.query_name, "New:" + soft_clip_seq, "Old:" + old_soft_clip_seq]
					break
			if flag == 0:
				out = ["FALSE", aln.query_name, "New:", "Old:" + old_soft_clip_seq]
		elif pattern == "TF":
			flag = 0
			for i in range( len(window_gc_list))[::-1]:
				if window_gc_list[i] == 0:
					flag = 1
					soft_clip_seq = soft_clip_seq[: (i + window_size) ]
					out = ["TRUE", aln.query_name, "New:" + soft_clip_seq, "Old:" + old_soft_clip_seq]
					break
			if flag == 0:
				out = ["FALSE", aln.query_name, "New:", "Old:" + old_soft_clip_seq]
		else:
			out = ["FAIL", aln.query_name, "New:" + soft_clip_seq, "Old:" + old_soft_clip_seq]

	# output
	out = '\t'.join([str(x) for x in out])
	print(out)
