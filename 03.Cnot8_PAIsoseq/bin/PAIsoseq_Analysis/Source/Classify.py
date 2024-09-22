#!/usr/bin/env python
import os
import pysam
from Source import Misc
import subprocess

class Classify():

	def __init__(self, parameters):

		print("\n-----------------")
		print("Classify CCS reads")
		print("------------------")

		# Output directory
		outputDir = parameters.getParametersBysubKey('experiment', 'outputDir')

		# Input file : Assignment result
		assignDir = os.path.join(outputDir, 'AssignDir')
		#self.assignPath = os.path.join(assignDir, 'clean_CCS.sorted.bam.featureCounts')
		self.assignPath = os.path.join(assignDir, 'clean_CCS.sorted.all_mapped.bam.featureCounts')
		Misc.Misc.checkFile(self.assignPath, 'Annotation of CCS reads')

		# Input file : input bad polyA file and alignment result
		aligntoTranscriptDir = os.path.join(outputDir, 'AligntoTranscriptDir')
		self.badPolya = os.path.join(aligntoTranscriptDir,'clean_CCS.sorted.all_mapped.badPolyA.refine.txt')
		Misc.Misc.checkFile(self.badPolya, 'badPolyA refine file')

		alignDir = os.path.join(outputDir, 'AligntoGenomeDir')
		self.bamPath = os.path.join(alignDir,'clean_CCS.sorted.all_mapped.bam')
		Misc.Misc.checkFile(self.bamPath, 'Alignment file')

		# Mitchondrial genome
		self.mtGenome = parameters.getParametersBysubKey('experiment', 'mtGenome')

		# chloroplast genome
		self.clGenome = parameters.getParametersBysubKey('experiment', 'clGenome')

        # Create output directory
		self.classifyDir = os.path.join(outputDir, 'ClassifyDir')
		Misc.Misc.checkDir(self.classifyDir, 'ClassifyDir')

	def classify(self):

		# Extract unmapped reads from bam file
#		print("\nExtract unmapped reads...")
#		unmappedReadPath = os.path.join(self.classifyDir, 'Unmapped.reads.txt')
#		print('Output: ' + unmappedReadPath)
#		unmappedReadFile = open(unmappedReadPath, 'w')
#		bamfile = pysam.AlignmentFile(self.bamPath, 'rb')
#		for aln in bamfile:
#			if aln.is_unmapped:
#				unmappedReadFile.write(aln.query_name + '\n')
#		unmappedReadFile.close()
		reads_list = []
		with open(self.badPolya, 'r') as badreads:
			for line in badreads:
				ID = line.strip().split("\t")[0]
				reads_list.append(ID)
		badreads.close()

		# MT reads
		# get mitchondrial chromosome name
		mtReadPath = os.path.join(self.classifyDir, 'MT.reads.txt')
		if self.mtGenome != "" and os.path.exists(self.mtGenome):
			print('\nExtract mitchondrial reads...')
			print('Output: ' + mtReadPath)
			mt_name = ""
			with open(self.mtGenome, 'r') as mt:
				for line in mt:
					line = line.strip('\n').split()
					mt_name = line[0][1:]
					break
			print("Mitchondrial genome name: " + str(mt_name))

			mtReadFile = open(mtReadPath, 'w')
			bamfile = pysam.AlignmentFile(self.bamPath, 'rb')
			for aln in bamfile:
				if aln.reference_name == mt_name and aln.query_name in reads_list:
					mtReadFile.write(aln.query_name + '\n')
			mtReadFile.close()

		# CL reads
		# get chloroplast chromosome name
		clReadPath = os.path.join(self.classifyDir, 'CL.reads.txt')
		if self.clGenome != "" and os.path.exists(self.clGenome):
			print('\nExtract chloroplast reads...')
			print('Output: ' + clReadPath)
			cl_name = ""
			with open(self.clGenome, 'r') as cl:
				for line in cl:
					line = line.strip('\n').split()
					cl_name = line[0][1:]
					break
			print("Chloroplast genome name: " + str(cl_name))

			clReadFile = open(clReadPath, 'w')
			bamfile = pysam.AlignmentFile(self.bamPath, 'rb')
			for aln in bamfile:
				if aln.reference_name == cl_name and aln.query_name in reads_list:
					clReadFile.write(aln.query_name + '\n')
			clReadFile.close()

		# Reads to exclude
		excludeReads = {}
		if self.mtGenome != "" and os.path.exists(mtReadPath):
			with open(mtReadPath, 'r') as mt:
				for line in mt:
					line = line.strip('\n').split()
					excludeReads[ line[0] ] = 1

		if self.clGenome != "" and os.path.exists(clReadPath):
			with open(clReadPath, 'r') as cl:
				for line in cl:
					line = line.strip('\n').split()
					excludeReads[ line[0] ] = 1

		# Assigned chromosome reads
		print("\nExtract assigned chromosome reads")
		AssignedReadsPath = os.path.join(self.classifyDir, 'Assigned.reads.txt')
		print('Output: ' + AssignedReadsPath)
		AssignedReadsFile = open(AssignedReadsPath, 'w')
		with open(self.assignPath, 'r') as fc:
			for line in fc:
				(read_name, status, strand, gene) = line.strip('\n').split()
				#if status == "Assigned" and excludeReads.has_key(read_name) == False:
				if status == "Assigned" and read_name in reads_list \
					and read_name not in excludeReads.keys():
					AssignedReadsFile.write(read_name + '\t' + gene + '\n')
		AssignedReadsFile.close()

		# Not assigned chromosome reads
		print("\nExtract not assigned chromosome reads")
		NotAssignedReadsPath = os.path.join(self.classifyDir, 'NotAssigned.reads.txt')
		print('Output: ' + NotAssignedReadsPath)
		NotAssignedReadsFile = open(NotAssignedReadsPath, 'w')
		with open(self.assignPath, 'r') as fc:
			for line in fc:
				(read_name, status, strand, gene) = line.strip('\n').split()
				if status != "Unassigned_Unmapped" and gene == "NA" \
					and read_name not in excludeReads.keys() \
					and read_name in reads_list:
					#and excludeReads.has_key(read_name) == False:
					NotAssignedReadsFile.write(read_name + '\n')
		NotAssignedReadsFile.close()
