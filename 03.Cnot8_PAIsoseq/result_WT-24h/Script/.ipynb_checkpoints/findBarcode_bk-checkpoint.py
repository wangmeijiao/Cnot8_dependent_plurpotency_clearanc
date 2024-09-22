#!/usr/bin/env python
import re
import sys
import pysam
import regex
import parasail
from itertools import islice

def revComp(s):
	revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
	return ''.join([revCDict[x] for x in s])[::-1]

def parseAlignment(result):
	# Traceback pair-wise alignment
	aligned_query  = result.traceback.query
	aligned_target = result.traceback.ref

	# Define gaps and alignment blocks
	regex = r"^(-*)([ATCG-]+?)(-*)$"
	matches = re.finditer(regex, aligned_query)
	(left_gap, mid, right_gap) = (0,0,0)
	for it in matches:
		left_gap = len(it.group(1))
		mid  = len(it.group(2))
		right_gap = len(it.group(3))
		break

	# Examples:
	# Query  : --------------ATCGAA-TAAA--------------
	# Target : AAAAAACTGAAAGCAT-GAATTAAA--------------
	# Group  : | Group1      |  Group2 | Group3      |
	# Notice : Group1 and Group3 may be empty

	# alignment start and end in reference/target
	aligned_target_nogap = aligned_target.replace('-','') # reference/target length
	target_len = len(aligned_target_nogap) # target sequence length
	start = left_gap # alignment start position in target
	end = target_len - right_gap # alignment end position in target

	# block coordinates
	aligned_query_block = aligned_query[left_gap:(left_gap + mid)]
	aligned_ref_block   = aligned_target[left_gap:(left_gap + mid)]

	# count match, mismatch, insert, deletion
	(match, mismatch, insert, deletion) = (0,0,0,0)
	for index in range(len(aligned_ref_block)):
		query_base = aligned_query_block[index]
		ref_base = aligned_ref_block[index]
		if query_base == ref_base: # Match
			match = match + 1
		elif query_base == "-":    # Insertion in reference
			insert = insert + 1
		elif ref_base == "-":      # Deletion in reference
			deletion = deletion + 1
		else:                      # Mismatch
			mismatch = mismatch + 1
	return match, mismatch, insert, deletion, start, end

def DetectAdapter(adapter, seq, qual, maxErr):
	'''
	Detect 5' TSO
	'''
	gap_open, gap_extend = 8, 4
	score_matrix = parasail.matrix_create("ACGT", 5, -5)
	result = parasail.sg_dx_trace(adapter, seq, gap_open, gap_extend, score_matrix)
	match, mismatch, insert, deletion, start, end = parseAlignment(result)
	error = mismatch + insert + deletion
	match_seq = ''
	match_qual = ''
	if error <= maxErr:
		match_seq = seq[start:end]
		match_qual = qual[start:end]
		seq = seq[end:]
		qual = qual[end:]
		#aligned_query  = result.traceback.query #add
		#aligned_target = result.traceback.ref # add
		#print("TSO     : " + aligned_query) # add
		#print("CCS Read: " + aligned_target) # add
	return (seq, qual, match_seq, match_qual)

def trimAdapter(adapter, seq, qual, maxErr):
	'''
	Remove adapters
	'''
	gap_open, gap_extend = 8, 4
	score_matrix = parasail.matrix_create("ACGT", 5, -5)
	while(True):
		result = parasail.sg_dx_trace(adapter, seq, gap_open, gap_extend, score_matrix)
		match, mismatch, insert, deletion, start, end = parseAlignment(result)
		error = mismatch + insert + deletion
		if seq != "" and error <= maxErr:
		#if error <= maxErr:
			seq = seq[end:]
			qual = qual[end:]
		else:
			break
	return (seq, qual)

def removeAdapter(seq, qual):
	'''
	Remove all possible TSO
	'''
	adapter = 'AAGCAGTGGTATCAACGCAGAGTACATGGG' # 25 nt, TSO
	(clean_cDNA, clean_qual) = trimAdapter(adapter, seq, qual, maxErr = 4)

	adapter = 'CCCATGTACTCTGCGTTGATACCACTGCTT'
	(clean_cDNA, clean_qual) = trimAdapter(adapter, clean_cDNA, clean_qual, maxErr = 4)

	adapter = 'AAGCAGTGGTATCAACGCAGAGTAC' # 25 nt, TSO
	(clean_cDNA, clean_qual) = trimAdapter(adapter, clean_cDNA, clean_qual, maxErr = 3)

	adapter = 'GTACTCTGCGTTGATACCACTGCTT' # 25 nt, TSO RC
	(clean_cDNA, clean_qual) = trimAdapter(adapter, clean_cDNA, clean_qual, maxErr = 3)

	adapter = 'TGCGTTGATACCACTGCTT' # 19 nt, TSO RC
	(clean_cDNA, clean_qual) = trimAdapter(adapter, clean_cDNA, clean_qual, maxErr = 2)

	adapter = 'AAGCAGTGGTATCAACGCA' # 19 nt, TSO RC
	(clean_cDNA, clean_qual) = trimAdapter(adapter, clean_cDNA, clean_qual, maxErr = 2)

	return (clean_cDNA, clean_qual)

def trimBarcode(barcode, seq, qual):
	'''
	Remove barcodes
	'''
	gap_open, gap_extend = 8, 4
	maxErr = 4
	barcode = revComp(barcode)
	seq = revComp(seq)
	score_matrix = parasail.matrix_create("ACGT", 5, -5)
	result = parasail.sg_dx_trace(barcode, seq,gap_open, gap_extend, score_matrix)
	#aligned_query  = result.traceback.query #add
	#aligned_target = result.traceback.ref # add
	#print("TSO     : " + revComp(aligned_query)) # add
	#print("CCS Read: " + revComp(aligned_target)) # add
	match, mismatch, insert, deletion, start, end = parseAlignment(result)
	error = mismatch + insert + deletion
	if seq != "" and error <= maxErr:
		seq = seq[end:]
		qual = qual[end:]
	seq = revComp(seq)
	qual = qual[::-1]
	return (seq, qual)


def cleanCCS(target_seq, target_qual, barcode_name, barcode_seq, mismatch, strand):
	'''
	Clean CCS Reads
	'''
	pattern = r"("+ barcode_seq + "){e<=" + str(mismatch) + "}"
	count = 0
	while(True):
		matches = regex.finditer( pattern, target_seq, overlapped = False)
		for it in matches:
			match_start = it.start()
			match_end = it.end()
			seq = target_seq[:match_end]
			qual = target_qual[:match_end]

			# Detect P5
			P5 = 'NO'
			p5Adapter = 'AAGCAGTGGTATCAACGCAGAGTACATGGG' 
                        #AAGCAGTGGTATCAACGCAGAGTAC
			(seq, qual, match_seq, match_qual) = DetectAdapter(p5Adapter, seq, qual, maxErr = 4)
			if match_qual != '':
				P5 = "YES"
			extract_seq = match_seq + seq
			extract_qual = match_qual + qual

			# Trim Adapter (trim all possible TSO sequence)
			(clean_seq, clean_qual) = removeAdapter(seq, qual)

			# Trim Barcode
			(clean_seq, clean_qual) = trimBarcode(barcode_seq[:16], clean_seq, clean_qual)

			# Output
			if clean_seq != "":
				count += 1
				read_pass = pass_dict[read.name]
				read_id = read.name + strand + str(count) + ':' + barcode_name + ':' + read_pass
				out = [read.name, read_id, barcode_name, P5 ,strand, extract_seq, extract_qual, clean_seq, clean_qual]
				out = '\t'.join([str(x) for x in out])
				if clean_seq[-8:].count('A') >= 4:
					print(out)
			else:
				sys.stderr.write(read.name + '\n')

			target_seq = target_seq[match_end:]
			target_qual = target_qual[match_end:]
			break
		else:
			break

if len(sys.argv) != 5:
	msg = ['python', sys.argv[0], 'infile', 'barcodefile','passfile','mismatch', '1>out.txt', '2>err.txt']
	msg = ' '.join([str(x) for x in msg])
	print(msg)
	print('\tinfile: Input raw CCS reads in fastq format')
	print('\tpassfile: CCS read pass file')
	print('\tbarcodefile: Barcode file in fasta format')
	print('\tmismatch: Number of mismatches')
	quit()

# Arguments
infile = sys.argv[1]
passfile = sys.argv[2]
barcodefile = sys.argv[3]
mismatch = int(sys.argv[4])

# Barcodes
barcode_dict = {}
with pysam.FastxFile(barcodefile) as fa_in:
	for cnt, read in enumerate(islice(fa_in,None)):
		barcode_dict[read.name] = read.sequence[:22] # 16 nt Barcode + 6 nt 3'-Adapter

# Pass file
pass_dict = {}
with open(passfile,'r') as pass_file:
	for line in pass_file:
		(read_name, read_pass) = line.strip('\n').split('\t')
		pass_dict[read_name] = read_pass

# Process
cts = 0
with pysam.FastxFile(infile) as fq_in:
	for cnt, read in enumerate(islice(fq_in,None)):
		cts += 1
		for barcode_name, barcode_seq in barcode_dict.items():
			target_seq = read.sequence
			target_qual = read.quality
			cleanCCS(target_seq, target_qual, barcode_name, barcode_seq, mismatch, strand = '+')

			target_seq = revComp(read.sequence)
			target_qual = read.quality[::-1]
			cleanCCS(target_seq, target_qual, barcode_name, barcode_seq, mismatch, strand = '-')

logfile = infile + '.log.txt'
log_file = open(logfile, 'w')
log_file.write(str(cts) + '\n')
log_file.close()
