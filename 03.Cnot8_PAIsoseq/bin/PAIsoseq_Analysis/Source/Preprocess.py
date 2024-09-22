#!/usr/bin/env python
import os
import sys
import regex
import pysam
import subprocess
from itertools import islice

class Preprocess():

	# Define Parameters for
	def __init__(self, parameters):

		# Get parameters
		self.parameters = parameters

		# Input files: CCS reads (fastq)
		self.fastqPath = parameters.getParametersBysubKey('experiment', 'rawFastq')

		# Input files: PASS file (txt)
		self.passFile = parameters.getParametersBysubKey('experiment', 'passFile')

		# Input files: Barcode file (fasta)
		self.barcodeFile = parameters.getParametersBysubKey('experiment', 'barcodeFile')

		# Check input Fastq File
		Misc.Misc.checkFile(self.fastqPath, 'Input Fastq File')
		Misc.Misc.checkFile(self.passFile, 'Pass file')
		Misc.Misc.checkFile(self.barcodeFile, 'Barcode file')

		# Check if outdir is registered, if so check PreprocessDir
		self.expOutdir = self.parameters.getParametersBysubKey('experiment', 'outputDir')

		# Generate preprocessing dir
		self.preprocessDir = os.path.join(self.expOutdir, 'PreprocessDir')
		Misc.Misc.checkDir(self.preprocessDir, 'PreprocessDir')

	def process(self):

		print('\n----------------')
		print("Preprocess Reads")
		print("----------------")

		# link fastq data
		fastqPath = os.path.join(self.preprocessDir, os.path.basename(self.fastqPath))
		if os.path.exists(fastqPath) == False:
			cmd = ['ln -s', self.fastqPath, self.preprocessDir]
			cmd = ' '.join([str(x) for x in cmd])
			print('\nLink reads (fastq):')
			print(cmd)
			os.system(cmd)
		self.fastqPath = os.path.abspath(fastqPath)

		# unzip
		if self.fastqPath.endswith('.gz'):
			cmd = ['gzip', '-dc', self.fastqPath, ">", self.fastqPath[ : -3]]
			cmd = ' '.join([str(x) for x in cmd])
			print('\nUnzip reads (fastq):')
			print(cmd)
			os.system(cmd)
			self.fastqPath = self.fastqPath[ : -3]

		# split fastq file into 40 parts
		cmd = ['fastq-splitter.pl', '--n-parts', 40, self.fastqPath]
		cmd = ' '.join([str(x) for x in cmd])
		print('\nSplit reads (fastq):')
		print(cmd)
		os.system(cmd)

		# collect sub files
		prefix = self.fastqPath.split('.')[:-1]
		prefix = '.'.join(prefix)
		suffix = self.fastqPath.split('.')[-1]
		parts = []
		for i in range(1,41):
			num = "{:0>2d}".format(i)
			parts.append(prefix + '.part-' + num + '.' + suffix)

		# demultiplex and trim barcode and adapter
		mismatch = 3
		run_file_path = os.path.join(self.preprocessDir, 'run.sh')
		run_file= open(run_file_path, 'w')
		for i in parts:
			out = '1>' + i + '.out'
			err = '2>' + i + '.err'
			cmd = ['findBarcode.py', i, self.passFile, self.barcodeFile, mismatch, out, err]
			cmd = ' '.join( [ str(x) for x in cmd])
			run_file.write(cmd + '\n')
		run_file.close()

		# run
		cmd = ['ParaFly', '-c', run_file_path, '-CPU', 10]
		cmd = ' '.join( [ str(x) for x in cmd])
		os.system(cmd)

		# collect results: out files
		out_files = []
		for part in parts:
			out_files.append( part + '.out')
		out_file = os.path.join(self.preprocessDir, 'clean_CCS.out.txt')
		cmd = ['cat']
		cmd.extend(out_files)
		cmd.extend(['1>' + out_file])
		cmd = ' '.join( [ str(x) for x in cmd])
		os.system(cmd)

		# collect results: err files
		err_files = []
		for part in parts:
			err_files.append( part + '.err')
		err_file = os.path.join(self.preprocessDir, 'clean_CCS.err.txt')
		cmd = ['cat']
		cmd.extend(err_files)
		cmd.extend(['1>' + err_file])
		cmd = ' '.join( [ str(x) for x in cmd])
		os.system(cmd)

		# collect results: log files
		log_files = []
		for part in parts:
			log_files.append(part + '.log.txt')
		log_file = os.path.join(self.preprocessDir, 'clean_CCS.log.txt')
		cmd = ['cat']
		cmd.extend(log_files)
		cmd.extend(['1>' + log_file])
		cmd = ' '.join( [ str(x) for x in cmd])
		os.system(cmd)

		# Clean temporary files
		cmd = ['rm']
		cmd.extend(parts)
		cmd.extend(out_files)
		cmd.extend(err_files)
		cmd.extend(log_files)
		cmd = ' '.join( [ str(x) for x in cmd])
		os.system(cmd)

		# Convert result to fasta and fastq
		fa_out_file = os.path.join(self.preprocessDir, 'clean_CCS.fa')
		fq_out_file = os.path.join(self.preprocessDir, 'clean_CCS.fq')
		fa_out = open(fa_out_file, 'w')
		fq_out = open(fq_out_file, 'w')
		with open(out_file, 'r') as ccsout:
			for line in ccsout:
				fields = line.strip('\n').split('\t')
				fa_record = [ ">" + fields[1], fields[-2] ]
				fa_record = '\n'.join([ str(x) for x in fa_record])
				fa_out.write(fa_record + '\n')

				fq_record = ["@" + fields[1], fields[-2], '+', fields[-1]]
				fq_record = '\n'.join( [str(x) for x in fq_record] )
				fq_out.write(fq_record + '\n')
		print("\nOutput fasta file:" + fa_out_file)
		print("\nOutput fastq file:" + fq_out_file)
		print("\nDone")
		ccsout.close()
		fa_out.close()
		fq_out.close()
