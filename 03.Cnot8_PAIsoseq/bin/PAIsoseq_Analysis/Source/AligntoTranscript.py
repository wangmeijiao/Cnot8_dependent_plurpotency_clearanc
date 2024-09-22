#!/usr/bin/env python3
import os
import subprocess
import pysam
from Source import Misc
import regex as re
from itertools import islice

class AligntoTranscript():

	def __init__(self, parameters):

		print("\n--------------------------")
		print("Align reads to transcripts")
		print("--------------------------")

		# Output directory
		outputDir = parameters.getParametersBysubKey('experiment', 'outputDir')

		# Input file
		alignGenomeDir = os.path.join(outputDir, 'AligntoGenomeDir')
		self.badpolyaPath = os.path.join(alignGenomeDir, \
						'clean_CCS.sorted.all_mapped.badPolyA.txt')
		Misc.Misc.checkFile(self.badpolyaPath, 'bad polyA file')

		# FeatureCounts file
		#assignDir = os.path.join(outputDir, 'AssignDir')
		#self.assignedPath = os.path.join(assignDir, \
		#				'clean_CCS.sorted.all_mapped.bam.featureCounts')
		#Misc.Misc.checkFile(self.assignedPath, 'featureCounts file')

        # Generate new dir for writing output files
		self.AlignTDir = os.path.join(outputDir, 'AligntoTranscriptDir')
		Misc.Misc.checkDir(self.AlignTDir, 'AligntoTranscriptDir')

		# Output files and log file
		self.badSCPath = os.path.join(self.AlignTDir, 'badPolya_3sc.txt')
		self.badSCLogPath = os.path.join(self.AlignTDir, 'badPolya_3sc.err')
		self.badSCWindowPath = os.path.join(self.AlignTDir, 'badPolya_3sc_window.txt')
		self.badSCFilterPath = os.path.join(self.AlignTDir, 'badPolya_3sc_window0.txt')
		self.badSCRefinePath = os.path.join(self.AlignTDir, \
						'clean_CCS.sorted.all_mapped.badPolyA.refine.txt')

        # Check if software is accressible
		# bwa
		self.bwaPath = parameters.getParametersBysubKey('software', 'bwa')

		# samtools
		self.samtoolsPath = parameters.getParametersBysubKey('software', 'samtools')

		# sam2pairwise
		self.sam2pairwisePath = parameters.getParametersBysubKey('software', 'sam2pairwise')

		if self.bwaPath != "":
			Misc.Misc.checkProg(self.bwaPath, 'bwa')

		Misc.Misc.checkProg(self.samtoolsPath, 'samtools')
		Misc.Misc.checkProg(self.sam2pairwisePath, 'sam2pairwise')

		# Number of threads used by the mapper : bwa
		self.nThreads = parameters.getParametersBysubKey('software', 'nThreads')

        # Check index (file or directory)
		self.cDNAindexPath = parameters.getParametersBysubKey('experiment', 'cDNAindex')
		Misc.Misc.checkFile(self.cDNAindexPath + '.ann', 'Check bwa index')
		Misc.Misc.checkFile(self.cDNAindexPath + '.amb', 'Check bwa index')
		Misc.Misc.checkFile(self.cDNAindexPath + '.bwt', 'Check bwa index')
		Misc.Misc.checkFile(self.cDNAindexPath + '.pac', 'Check bwa index')
		Misc.Misc.checkFile(self.cDNAindexPath + '.sa', 'Check bwa index')

		self.thread = 10
		self.windowsize = 10
		self.stepsize = 5
		self.extendPolyA = False

	def aligntoTranscript(self):
		"""
		Align bad polyA CCS reads to reference cDNA using bwa and 
		get the confident poly(A) tails.
		"""

		# Get bad polyA reads
		badReadPath = os.path.abspath(os.path.join(self.AlignTDir, 'badpolyA_reads.fa'))
		cmd = ['cat', self.badpolyaPath, "| awk '{print \">\"$1\"\\n\"$NF}' >",badReadPath]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nExtract bad polyA reads:")
		print(cmd)
		os.system(cmd)

		#Mapping
		out_bam = os.path.abspath(os.path.join(self.AlignTDir,'badPolyA_reads.bam'))
		err_txt = os.path.abspath(os.path.join(self.AlignTDir, 'align.log'))
		out = '| ' + self.samtoolsPath + ' view -bS -F3840 | samtools sort - >' + out_bam
		err = '2> ' + err_txt
		# align cmd
		#-A2 -B1 -O1 -E1
		cmd = [self.bwaPath, "mem", self.cDNAindexPath, badReadPath, '-t', self.nThreads, \
		'-k11 -W40 -r10 -A2 -B1 -O4 -E2 -L0 -T 20', err, out ]
		#cmd = [self.bwaPath, "mem", self.cDNAindexPath, badReadPath, '-t', self.nThreads, err, out ]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nMap bad polyA reads to cDNA:")
		print(cmd)
		os.system(cmd)

		badSC = open(self.badSCPath, 'w')
		badSCFilter = open(self.badSCFilterPath, 'w')

		# extract extended 3'soft-clip
		dic_sc,refine_list = self.extract_soft_clip(out_bam)

		for k in dic_sc.keys():
			three_soft_clip_seq = dic_sc[k]
			read_ID = k

			#filter G/C in the start of three_soft_clip_seq
			#if three_soft_clip_seq != "":
			#	while(three_soft_clip_seq[0] == "C" or three_soft_clip_seq[0] == "G"):
			#		three_soft_clip_seq = three_soft_clip_seq[1:]
			#		if three_soft_clip_seq == "":
			#			break

			if three_soft_clip_seq == '':
				continue
			else:
				badSC.write(read_ID + "\t" + three_soft_clip_seq+ "\n")		
			# Whether to extend A
			# Extend A when the first residue of the 3'-soft clip is A

			if three_soft_clip_seq != "":
			#	while( 'C' in three_soft_clip_seq[:3] or 'G' in three_soft_clip_seq[:3]):
			#		three_soft_clip_seq = three_soft_clip_seq[1:]
			#	if three_soft_clip_seq != "":
					# filter bad polyA by window score
				#window = self.windowsize
				#step_size = self.stepsize
				ct = 0
				window_score = ""
				for i in range( 0, len(three_soft_clip_seq), self.stepsize):
					window_seq = three_soft_clip_seq[i : (i + self.windowsize) ]
					ct += 1
					if window_seq.count('A')  < int(0.7*len(window_seq)):
						window_score = window_score + "F"
					else:
						window_score = window_score + "T"
					if len(window_seq) < self.windowsize:
						break
				#if window_score.count("F") <= 1 and \
				#	three_soft_clip_seq[:3].count('A') == len(three_soft_clip_seq[:3]):
				#	badSCFilter.write(aln.query_name + "\t" + three_soft_clip_seq+ "\n")	
				if window_score.count('F') == 0 :#and \
					#three_soft_clip_seq[:2].count('A') == len(three_soft_clip_seq[:2]):
					while "G" in three_soft_clip_seq[:2] or "C" in three_soft_clip_seq[:2]:
						three_soft_clip_seq = three_soft_clip_seq[1:]
						if three_soft_clip_seq == "":
							break
					#	three_soft_clip_seq = three_soft_clip_seq[1:]
					#	if three_soft_clip_seq == "":
					#		break
					if three_soft_clip_seq != "":
						badSCFilter.write(read_ID + "\t" + three_soft_clip_seq+ "\n")
					else:
						refine_list.append(read_ID)
				else:
					refine_list.append(read_ID)
		badSC.close()
		badSCFilter.close()

		# Get remaining reads for next step
		badSCRefine = open(self.badSCRefinePath, 'w')
		for ID in refine_list:
			badSCRefine.write(ID + "\n")
		badSCRefine.close()
		cmd = ['fisher.pl', self.badSCRefinePath, self.badpolyaPath, "> tmp && mv tmp", self.badSCRefinePath]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nExtract remaining bad polyA reads:")
		print(cmd)
		os.system(cmd)

	def extract_soft_clip(self,bam):
		prefix = bam.rstrip('.bam')
		pairwiseFile =  prefix + '.sam2pairwise.out'
		print("\nConvert bam format to pairwise format:")
		cmd = [self.samtoolsPath,'view',bam, '|', self.sam2pairwisePath, '>',pairwiseFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
		os.system(cmd)

		dic_sc = {}
		refine_list = []
		tail = ""

		with open(pairwiseFile, 'r') as pwfile:
			while True:
				lines = list(islice(pwfile, 4))
				if not lines:
					break
				inform,read,char,ref = lines
				inform = inform.strip()
				read = read.strip()
				char = char.strip()
				ref = ref.strip()
				ID = inform.split("\t")[0]
				flag = inform.split("\t")[1]
				n=0
				if flag == "0":
					for i in ref[::-1]:
						if i == "N":
							n += 1
						else:
							break
					for j in range(-(n+1),-(len(ref)+1),-1):
						if abs(j) == len(ref):
							read = read.replace('-','')
							read = read.replace('.','')
							tail = read
						elif read[j] != "A" and ref[j] != "A" and read[j] == ref[j]:
							tail = read[j+1:]
							tail = tail.replace('-','')
							tail = tail.replace('.','')
							break
						else:
							continue
					if tail != "":
						dic_sc[ID] = tail
					else:
						refine_list.append(ID)
				else:
					refine_list.append(ID)
		return dic_sc,refine_list

	def revComp(self,s):
		revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
		return ''.join([revCDict[x] for x in s])[::-1]

