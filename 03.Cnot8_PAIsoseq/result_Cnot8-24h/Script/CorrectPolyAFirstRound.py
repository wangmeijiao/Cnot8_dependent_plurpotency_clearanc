#!/usr/bin/env python
import os
import sys
import pysam
from itertools import islice

if len(sys.argv) != 4:
	msg = ['python',sys.argv[0],'soft_clip_file','polyA_file', 'down', '1>out', '2>err']
	msg = ' '.join([str(x) for x in msg])
	print(msg)
	print("\tRound 1, Correct polyA tail")
	print('\tsoft_clip_file: Soft clip file (fasta)')
	print('\tpolyA_file: Draft polyA file')
	print('\tdown: 10')
	print('\toutfile: Stand out')
	quit()

infile = sys.argv[1]
polyAfile = sys.argv[2]
down = int(sys.argv[3])

# read soft clip
sc_dict ={}
with pysam.FastxFile(infile) as fa_file:
	for cnt, e in enumerate(islice(fa_file,None)):
		sc_dict[e.name] = e.sequence

# Record
keep = 0
correct = 0

logfile = polyAfile + '.correct.log'
with open(polyAfile,'r') as polyAfile:
	for line in polyAfile:
		flag = '' # Corrected: YES? or NO?
		edit_num = 0
		extendA_num = 0

		(name, polyA_seq) = line.strip('\n').split('\t')
		soft_clip = 0
		if name in sc_dict:
			soft_clip = len(sc_dict[name])
		else:
			sys.stderr.write(name+'\n')
			continue

		trim_num = 0
		if len(polyA_seq) <= down: # polyA length <= 10
			if len(polyA_seq) >= soft_clip:
				flag = 'YES'
				correct += 1
				trim_num = len(polyA_seq) - soft_clip
				edit_num = trim_num

				trimmed_seq = polyA_seq[trim_num:]
				for idx in range(trim_num)[::-1]:
					if polyA_seq[idx] == 'A':
						trimmed_seq = 'A' + trimmed_seq
					else:
						break
				polyA_seq = trimmed_seq

			else: # soft clip length > polyA length
				flag = 'NO'
				edit_num = 0###edit_num = soft_clip - len(polyA_seq)
				extendA_num = 0
				keep += 1
				# change***
				###polyA_seq = sc_dict[name]

		else: # polyA length > 10
			if down >= soft_clip: # soft clip length <= 10
				flag =  'YES'
				correct += 1
				trim_num = down - soft_clip
				edit_num = trim_num

				trimmed_seq = polyA_seq[trim_num:]
				for idx in range(trim_num)[::-1]:
					if polyA_seq[idx] == 'A':
						extendA_num += 1
						trimmed_seq = 'A' + trimmed_seq
					else:
						break
				polyA_seq = trimmed_seq

			else: # soft clip length >= 10
				flag = 'NO'
				edit_num = 0 #####edit_num = soft_clip - down
				extendA_num = 0
				keep += 1
				# change***
				#####polyA_seq = sc_dict[name] + polyA_seq[10:]

		out = [name, flag, edit_num, extendA_num, polyA_seq]
		out = '\t'.join([str(x) for x in out])
		print(out)

with open(logfile, 'w') as log:
	log.write('summary:\nCorrected: ' + str(correct) + '\nUncorrected: ' + str(keep))
