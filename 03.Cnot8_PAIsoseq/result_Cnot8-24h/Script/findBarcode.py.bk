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

def trimAdapter(adapter, seq, qual, maxErr):
	'''
	Remove adapters
	'''
	gap_open, gap_extend = 8, 4
	score_matrix = parasail.matrix_create("ACGT", 5, -5)
	while(True):
		if seq == "":
			break
		result = parasail.sg_dx_trace(adapter, seq, gap_open, gap_extend, score_matrix)
		match, mismatch, insert, deletion, start, end = parseAlignment(result)
		error = mismatch + insert + deletion
		if error <= maxErr:
			seq = seq[end:]
			qual = qual[end:]
		else:
			break
	return (seq, qual)

def removeAdapter(seq, qual):
	'''
	Remove all possible TSO
	'''
	adapter = 'AAGCAGTGGTATCAACGCAGAGTAC'
	(clean_cDNA, clean_qual) = trimAdapter(adapter, seq, qual, maxErr = 3)

	adapter = 'GTACTCTGCGTTGATACCACTGCTT'
	(clean_cDNA, clean_qual) = trimAdapter(adapter, clean_cDNA, clean_qual, maxErr = 3)

	return (clean_cDNA, clean_qual)

def trimBarcode(barcode, seq, qual):
	'''
	Remove barcodes
	'''
	gap_open, gap_extend = 8, 4
	maxErr = 2
	barcode = revComp(barcode)
	seq = revComp(seq)
	score_matrix = parasail.matrix_create("ACGT", 5, -5)
	result = parasail.sg_dx_trace(barcode, seq,gap_open, gap_extend, score_matrix)
	match, mismatch, insert, deletion, start, end = parseAlignment(result)
	error = mismatch + insert + deletion
	if seq != "" and error <= maxErr:
		seq = seq[end:]
		qual = qual[end:]
	seq = revComp(seq)
	qual = qual[::-1]
	return (seq, qual)

def cleanCCS(transcript, mismatch, count):
	'''
	Clean CCS Reads
	'''
	(barcode_name, barcode_seq, target_name, target_seq, target_qual, strand) = transcript

	UMI = ""
	pattern = r"("+ barcode_seq + "){e<=" + str(mismatch) + "}"

	# Detect P3 (barcode + adapter)
	matches = regex.finditer( pattern, target_seq, overlapped = False)
	for it in matches:
		match_start = it.start()
		match_end = it.end()
		seq = target_seq[:match_end]
		qual = target_qual[:match_end]

		# Detect P5 (adapter + UMI + ATGGG)
		p5Adapter = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=2}([ATCG]{9,11})(ATGGG){e<=1}'
		p5_matches = regex.finditer(p5Adapter, seq, overlapped = False)
		for p5_it in p5_matches:
			(p5_tso, UMI,atggg) = p5_it.groups()
			p5_match_start = p5_it.start()
			p5_match_end = p5_it.end()
			seq  = seq[p5_match_end:]
			qual = qual[p5_match_end:]
			break # end
		else: # No P5 detected
			err = [target_name, barcode_name, strand, 'NO_P5']
			err = '\t'.join([str(x) for x in err])
			sys.stderr.write(err + '\n')
			break

		# Trim barcode
		(seq, qual) = trimBarcode(barcode_seq[:16], seq, qual)

		# Trim adapter
		(clean_seq, clean_qual) = removeAdapter(seq, qual)

		# Output
		if(clean_seq != ""):
			read_pass = pass_dict[target_name]
			read_id = target_name + strand + str(count) + ':' + barcode_name + ':' + read_pass
			out = [target_name, read_id, barcode_name, UMI ,strand, clean_seq, clean_qual]
			out = '\t'.join([str(x) for x in out])
			print(out)
		else:
			sys.stderr.write( target_name + '\n')

		target_seq  = target_seq[match_end:]
		target_qual = target_qual[match_end:]
		break # end

def locateBarcode(target_name, target_seq, barcode_name, barcode_seq, mismatch):
	'''
	Locate barcode in CCS reads
	'''
	location = {}
	pattern = r"("+ barcode_seq + "){e<=" + str(mismatch) + "}" # # 16 nt Barcode + 6 nt 3'-Adapter
	matches = regex.finditer( pattern, target_seq, overlapped = False)
	for it in matches: # barcode name, match start, match end
		location[it.start()] = (barcode_name, it.start(), it.end())
	return location

def splitCCS(read, barcode_dict, mismatch, strand):
	'''
	Split ccs by barcode
	'''
	# find barcode location
	all_bc_location = {}
	for barcode_name, barcode_seq in barcode_dict.items():
		if strand == "+":
			target_seq = read.sequence
			target_qual = read.quality
		elif strand == "-":
			target_seq = revComp(read.sequence)
			target_qual = read.quality[::-1]

		bc_location = locateBarcode(read.name, target_seq, barcode_name, barcode_seq, mismatch)
		all_bc_location.update(bc_location)

	# sort barcode location
	all_bc_sorted_location = []
	for key in sorted(all_bc_location.keys()):
		all_bc_sorted_location.append( all_bc_location[key] )

	# split ccs reads by barcode
	splits = []
	last_end = 0
	for i in range(len(all_bc_sorted_location)):
		this_end = all_bc_sorted_location[i][2]
		bc = all_bc_sorted_location[i][0]
		bc_seq = barcode_dict[bc]

		seq = target_seq[last_end:this_end]
		qual = target_qual[last_end:this_end]

		transcript = (bc, bc_seq, read.name, seq, qual, strand)
		splits.append(transcript)
		last_end = this_end

	return splits

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
cts = 0 # record total number of CCS reads
with pysam.FastxFile(infile) as fq_in:
	for cnt, read in enumerate(islice(fq_in,None)):
		cts += 1
		# split ccs by barcode
		splits_plus  = splitCCS(read, barcode_dict, mismatch, "+") # plus strand
		splits_minus = splitCCS(read, barcode_dict, mismatch, "-") # minus strand
		splits = splits_plus
		splits.extend(splits_minus)

		if len(splits) == 0: # No barcode found in CCS read
			sys.stderr.write(read.name + '\tNO_Barcode\n')
			continue

		# extract transcript
		count = 0
		for transcript in splits:
			count = count + 1
			cleanCCS(transcript, mismatch, count)

logfile = infile + '.log.txt'
log_file = open(logfile, 'w')
log_file.write(str(cts) + '\n')
log_file.close()
