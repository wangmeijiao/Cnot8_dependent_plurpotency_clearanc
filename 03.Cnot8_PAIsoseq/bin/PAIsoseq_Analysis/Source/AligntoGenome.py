#!/usr/bin/env python3
import os
import subprocess
import pysam
import regex as re
from Source import Misc
from itertools import islice

class AligntoGenome():

	def __init__(self, parameters):

		print("\n--------------------------------")
		print("Align reads to reference genomes")
		print("--------------------------------")

		# Output directory
		outputDir = parameters.getParametersBysubKey('experiment', 'outputDir')

		# Input file
		preprocessDir = os.path.join(outputDir, 'PreprocessDir')
		self.cleanCCSFqPath = os.path.join(preprocessDir, 'clean_CCS.fa')
		Misc.Misc.checkFile(self.cleanCCSFqPath, 'Clean CCS reads in fastq format')

        # Generate new dir for writing output files
		self.AligntoGenomeDir = os.path.join(outputDir, 'AligntoGenomeDir')
		Misc.Misc.checkDir(self.AligntoGenomeDir, 'AligntoGenomeDir')

		# Output sam/bam files and log file
		self.outBamPath = os.path.join(self.AligntoGenomeDir, 'clean_CCS.bam')
		self.outSamPath = os.path.join(self.AligntoGenomeDir, 'clean_CCS.sam')
		self.sortedBamPath = os.path.join(self.AligntoGenomeDir, 'clean_CCS.sorted.bam')
		self.alignLog = os.path.join(self.AligntoGenomeDir, 'align.log')

        # Check if software is accressible
		# deSALT
		self.deSALTPath = parameters.getParametersBysubKey('software', 'deSALT')
		self.loadGtfPath = parameters.getParametersBysubKey('software', 'loadGTF')

		# minimap2
		self.minimap2Path = parameters.getParametersBysubKey('software', 'minimap2')
		self.paftoolsPath = parameters.getParametersBysubKey('software', 'paftools')

		# STARlong
		self.STARlongPath = parameters.getParametersBysubKey('software', 'STARlong')

		# samtools, bedtools, bwa
		self.samtoolsPath = parameters.getParametersBysubKey('software', 'samtools')
		self.bedtoolsPath = parameters.getParametersBysubKey('software', 'bedtools')
		self.bwaPath = parameters.getParametersBysubKey('software', 'bwa')
		self.sam2pairwisePath = parameters.getParametersBysubKey('software', 'sam2pairwise')

		if self.deSALTPath != "":
			Misc.Misc.checkProg( self.deSALTPath, 'deSALT')
		if self.loadGtfPath != "":
			Misc.Misc.checkProg(self.loadGtfPath, 'Annotation_Load.py')

		if self.minimap2Path != "":
			Misc.Misc.checkProg(self.minimap2Path, 'minimap2')
		if self.paftoolsPath != "":
			Misc.Misc.checkProg(self.paftoolsPath, 'paftools.js')

		if self.STARlongPath != "":
			Misc.Misc.checkProg(self.STARlongPath, 'STARlong')

		Misc.Misc.checkProg(self.samtoolsPath, 'samtools')
		Misc.Misc.checkProg(self.bedtoolsPath, 'bedtools')
		Misc.Misc.checkProg(self.bwaPath, 'bwa')
		Misc.Misc.checkProg(self.sam2pairwisePath, 'sam2pairwise')

		# Number of threads used by the mapper : deSALT/minimap2/STARlong
		self.nThreads = parameters.getParametersBysubKey('software', 'nThreads')

		# Check BWA index, five files
		self.bwaIndex = parameters.getParametersBysubKey('experiment', 'bwaIndex')
		Misc.Misc.checkFile(self.bwaIndex + '.ann', 'Check bwa index')
		Misc.Misc.checkFile(self.bwaIndex + '.amb', 'Check bwa index')
		Misc.Misc.checkFile(self.bwaIndex + '.bwt', 'Check bwa index')
		Misc.Misc.checkFile(self.bwaIndex + '.pac', 'Check bwa index')
		Misc.Misc.checkFile(self.bwaIndex + '.sa', 'Check bwa index')

        # Check index (file or directory)
		self.genomeIndexPath = parameters.getParametersBysubKey('experiment', 'genomeIndex')
		if not os.path.exists(self.genomeIndexPath):
			sys.stderr.write("Can not locate index at {}\n".format(self.genomeIndexPath))
			quit()

		# Check genome annotation file
		self.gtfPath = parameters.getParametersBysubKey('experiment', 'annotationGTF')
		if self.gtfPath != "":
			Misc.Misc.checkFile(self.gtfPath, 'Genome annotations')

		# whether to extend polyA
		self.extendPolyA = False

	def aligntoGenome(self):
		if self.deSALTPath!= "":
			# Run deSALT
			self.alignReads(mapper = "deSALT")

		elif self.minimap2Path != "":
			# Run minimap2
			self.alignReads(mapper = "minimap2")

		elif self.STARlongPath != "":
			# Run STARlong
			self.alignReads(mapper = "STARlong")

	def alignReads(self, mapper):
		"""
		Align clean CCS reads to reference genomes using deSALT/minimap2/STARlong
		with or without gene annotation, and get the confident poly(A) tails.
		"""
		# Generate commands for mappers
		cmd = []
		if mapper == "deSALT":      # use deSALT
			cmd1 = [] 				# load gtf file
			cmd2 = [] 				# align CCS reads to reference genome
			cmd3 = [] 				# convert sam to bam
			cmd4 = [] 				# delete sam
			if self.gtfPath == "":
				cmd2 = [self.deSALTPath, 'aln', '-N 1 -x ccs', '-t', self.nThreads, \
						'-d 10 -T -f first_pass', '-o', self.outBamPath, \
						self.genomeIndexPath, self.cleanCCSFqPath,'&>' + self.alignLog]
			else:
				genomeInfoPath = os.path.join( os.path.abspath(self.AligntoGenomeDir),'genome.info')
				cmd1 = [self.loadGtfPath, self.gtfPath, genomeInfoPath]
				cmd2 = [self.deSALTPath, 'aln', '-N 1 -x ccs', '-t', self.nThreads, \
				'-d 10 -T -f first_pass','-G', genomeInfoPath, '-o', self.outSamPath, \
				self.genomeIndexPath, self.cleanCCSFqPath, "&>", self.alignLog]

			# Convert SAM to BAM
			cmd3 = [ self.samtoolsPath, 'view -bS', self.outSamPath, '>', self.outBamPath]

			# Delete SAM file (just keep BAM file)
			cmd4 = [ 'rm -rf', self.outSamPath]

			# combine command
			if self.gtfPath == "":
				cmd.extend(cmd2)
				cmd.append('&&')
				cmd.extend(cmd3)
				cmd.append('&&')
				cmd.extend(cmd4)
			else:
				cmd.extend(cmd1)
				cmd.append('&&')
				cmd.extend(cmd2)
				cmd.append('&&')
				cmd.extend(cmd3)
				cmd.append('&&')
				cmd.extend(cmd4)
	
		elif mapper == "minimap2": # use Minimap2
			if self.gtfPath == "":
				cmd = [self.minimap2Path, '-ax splice -uf --secondary=no --MD --cs', \
				'-t', self.nThreads, self.genomeIndexPath, self.cleanCCSFqPath, \
				' | ' + self.samtoolsPath + ' view -bS > ', self.outBamPath, \
				'2>' + self.alignLog]
			else:
				juncBedPath = os.path.join( os.path.abspath(self.AligntoGenomeDir),'junctions.bed')
				cmd1 = [self.paftoolsPath, 'gff2bed', self.gtfPath, '>', juncBedPath]
				cmd2 = [self.minimap2Path, '-ax splice -uf --secondary=no --MD --cs', \
				'-t', self.nThreads, '--junc-bed', juncBedPath, self.genomeIndexPath, \
				self.cleanCCSFqPath, ' | ' + self.samtoolsPath + ' view -bS > ', \
				self.outBamPath, '2>' + self.alignLog]

				cmd.extend(cmd1)
				cmd.append('&&')
				cmd.extend(cmd2)

		elif mapper == "STARlong": # use STARlong
			outPrefix = os.path.join( os.path.abspath(self.AligntoGenomeDir),'tmp')
			cmd1 = [ self.STARlongPAth, '--runThreadN', self.nThreads, \
			'--genomeDir', self.genomeIndexPath, '--readFiles', self.cleanCCSFqPath, \
			'--outFileNamePrefix', outPrefix, '--outSAMtype BAM SortedByCoordinate',\
			'--outFilterMultimapScoreRange 20, --outFilterScoreMinOverLread 0',\
			'--outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 1000', \
			'--winAnchorMultimapNmax 200 --seedSearchStartLmax 12 --seedPerReadNmax 100000',
			'--seedPerWindowNmax 100 --alignTranscriptsPerReadNmax 100000',\
			'--alignTranscriptsPerWindowNmax 10000 --readNameSeparator space', \
			'&>', self.alignLog ]
			cmd2 = ['mv', outPrefix + 'Aligned.sortedByCoord.out.bam', self.outBamPath]
			cmd.extend(cmd1)
			cmd.append('&&')
			cmd.extend(cmd2)

		# Generate command to align CCS reads to reference genomes
		# with or without gene annotation
		cmd = ' '.join([ str(x) for x in cmd])
		print("\nAlign clean ccs reads to reference genome:")
		print(cmd)
#		os.system(cmd)

		# sort
		print("\nSort BAM file:")
		cmd = [ self.samtoolsPath, 'sort', '-o', self.sortedBamPath, self.outBamPath]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

		# index
		print("\nIndex Sorted BAM file:")
		cmd = [self.samtoolsPath, 'index', self.sortedBamPath]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

		# extract unmapped reads and mapped reads
		prefix = self.sortedBamPath.rstrip('.bam')
		unmappedBamFile = prefix + '.unmapped.bam'
		print("\nExtract unmapped reads:")
		cmd = [self.samtoolsPath, "view -bS -f4", self.sortedBamPath, '>', unmappedBamFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

		# convert bam to fastq (unmapped reads)
		unmappedFqFile = prefix + '.unmapped.fq'
		cmd = [self.bedtoolsPath, 'bamtofastq', '-i', unmappedBamFile,'-fq', unmappedFqFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

		# extract mapped bam
		mappedBamFile = prefix + '.mapped.bam'
		print("\nExtract mapped BAM:")
		cmd = [self.samtoolsPath, 'view -bS -F3844', self.sortedBamPath, '>', mappedBamFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

		# align to reference using BWA MEM
		bwaBamFile = prefix + '.bwa_mapped.bam'
		print("\nAlign unmapped reads to reference genome using BWA MEM:")
		cmd = [self.bwaPath, 'mem', self.bwaIndex, unmappedFqFile, '-t', self.nThreads, \
		'-k17 -W40 -r10 -A1 -B2 -O1 -E1 -L0', '-T 20 | ', self.samtoolsPath, \
		'view -bS -F3844 | ',self.samtoolsPath,'sort - >', bwaBamFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

		# combine deSALT/minimap2/STARlong mapped reads and BWA mapped reads
		allMappedBamFile = prefix + '.all_mapped.bam'
		print("\nCombined two mapped BAM files:")
		cmd = [self.samtoolsPath,'merge -f', allMappedBamFile, mappedBamFile, bwaBamFile, 
		      '&&', self.samtoolsPath, 'index', allMappedBamFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

		#extend polyA tails
		print("\nPick nice polyA tails...")
		self.pick_nice_polyA(allMappedBamFile)
		print("Done")

	def pick_nice_polyA(self,bam):
		prefix = bam.rstrip('.bam')
		pairwiseFile =  prefix + '.sam2pairwise.out'
		print("\nConvert bam format to pairwise format:")
		cmd = [self.samtoolsPath,'view',bam, '|', self.sam2pairwisePath, '>',pairwiseFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)
		
		noPolyaFile = prefix + '.noPolyA.txt'
		nicePolyaFile = prefix + '.nicePolyA.txt'
		badPolyaFile = prefix + '.badPolyA.txt'

		noPolya = open(noPolyaFile, 'w')
		nicePolya = open(nicePolyaFile, 'w')
		badPolya = open(badPolyaFile, 'w')

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
						if read[j] != "A" and ref[j] != "A" and read[j] == ref[j]:
							tail = read[j+1:]
							tail = tail.replace('-','')
							tail = tail.replace('.','')
							read = read.replace('-','')
							read = read.replace('.','')
							break
						elif abs(j) == len(ref):
							read = read.replace('-','')
							read = read.replace('.','')
							tail = read
						else:
							continue
				else:
					for i in ref:
						if i == "N":
							n += 1
						else:
							break
					for j in range(n,len(ref),1):
						if read[j] != "T" and ref[j] != "T" and read[j] == ref[j]:
							tail = read[:j]
							tail = tail.replace('-','')
							tail = tail.replace('.','')
							tail = self.revComp(tail)
							read = read.replace('-','')
							read = read.replace('.','')
							read = self.revComp(read)
							break
						elif j == len(ref) - 1:
							read = read.replace('-','')
							read = read.replace('.','')
							read = self.revComp(read)
							tail = read
						else:
							continue


				# no polyA tail
				# The CCS read has no polyA tail:
				# 1) the soft clip length is 0,
				# 2) the soft clip length > 0, but not contain A residue,
				if 'A' not in tail:
					noPolya.write(ID + "\n")
				else:
					# define nice polyA
					if 'T' not in tail and \
					  'G' not in tail and \
				 	 'C' not in tail:
						nicePolya.write(ID+'\t'+tail+'\n')
					else:
						updist = 35
						badPolya.write(ID + '\t' + tail + '\t' + read[ (len(read) - len(tail) - updist): ] + "\t" + read +'\n')
		nicePolya.close()
		badPolya.close()
		noPolya.close()

		cmd = ['rm',pairwiseFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
#		os.system(cmd)

	def revComp(self,s):
		revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
		return ''.join([revCDict[x] for x in s])[::-1]
