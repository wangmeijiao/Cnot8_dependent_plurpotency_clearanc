#!/usr/bin/env python
import os
from Source import Misc
import subprocess
import pysam
from itertools import islice

class Refine():

	def __init__(self, parameters):

		print("\n------------------")
		print("Refine poly(A) tail")
		print("------------------")

		# Output directory
		outputDir = parameters.getParametersBysubKey('experiment', 'outputDir')
		outputDir = os.path.abspath(outputDir)

		# Chromosoeme genomes
		self.genomeFasta = parameters.getParametersBysubKey('experiment', 'genomeFasta')
		Misc.Misc.checkFile(self.genomeFasta, 'Chromosome genome fasta')

		# Gene index
		self.geneIndex = parameters.getParametersBysubKey('experiment', 'geneIndex')
		Misc.Misc.checkDir(self.geneIndex, 'Gene index directory')

		# Input file: draft polyA tail
		draftPolyaDir = os.path.join(outputDir, 'DraftPolyaDir')
		self.polyaPath= os.path.abspath(os.path.join(draftPolyaDir, 'clean_CCS_sorted_all_mapped_badPolyA_refine_draft.txt'))
		Misc.Misc.checkFile(self.polyaPath, 'PolyA tail information')

		# Input files: classification of CCS reads
		classifyDir = os.path.join(outputDir, 'ClassifyDir')
		self.mtReads = os.path.abspath(os.path.join(classifyDir, 'MT.reads.txt'))
		self.clReads = os.path.abspath(os.path.join(classifyDir, 'CL.reads.txt'))
		self.assignedReads = os.path.abspath(os.path.join(classifyDir, 'Assigned.reads.txt'))
		self.notassignedReads = os.path.abspath(os.path.join(classifyDir, 'NotAssigned.reads.txt'))

		# Input file: alignment bam (To deal with not assigned reads)
		alignDir = os.path.join(outputDir, 'AligntoGenomeDir')
		self.bamFile = os.path.abspath(os.path.join(alignDir, 'clean_CCS.sorted.all_mapped.bam'))

		# Input file: Clean CCS reads (cDNA include polyA tail)
		self.cleanCcsPath = os.path.join(outputDir, 'PreprocessDir', 'clean_CCS.fa')

        # Generate new dir for writing output files
		self.refinePolyaDir = os.path.abspath(os.path.join(outputDir, 'RefinePolyaDir'))
		Misc.Misc.checkDir(self.refinePolyaDir, 'RefinePolyaDir')

		# Mitchondrial genome
		self.mtGenome = parameters.getParametersBysubKey('experiment', 'mtGenome')

		# Chloroplast genome
		self.clGenome = parameters.getParametersBysubKey('experiment', 'clGenome')

		# Assigned directory
		self.assignedDir = os.path.abspath(os.path.join(self.refinePolyaDir, 'Assigned'))

		# Not Assigned directory
		self.notAssignedDir= os.path.abspath(os.path.join(self.refinePolyaDir, 'NotAssigned'))

        # Check if software is accressible
		self.bedtoolsPath = parameters.getParametersBysubKey('software', 'bedtools')
		Misc.Misc.checkProg(self.bedtoolsPath, 'bedtools')

		self.bwaPath = parameters.getParametersBysubKey('software', 'bwa')
		Misc.Misc.checkProg(self.bwaPath, 'bwa')

		self.samtoolsPath = parameters.getParametersBysubKey('software', 'samtools')
		Misc.Misc.checkProg(self.samtoolsPath, 'samtools')

		self.seqkitPath = parameters.getParametersBysubKey('software', 'seqkit')
		Misc.Misc.checkProg(self.seqkitPath, 'seqkit')

		self.sam2pairwisePath = parameters.getParametersBysubKey('software', 'sam2pairwise')
		Misc.Misc.checkProg(self.sam2pairwisePath, 'sam2pairwise')

		self.thread = 10

	def refine(self):

		# Correct MT reads
		if self.mtGenome != "" and os.path.exists(self.mtGenome) and os.path.getsize(self.mtReads) >0:
			self.refine_organelle("MT", self.mtGenome, self.mtReads)

		# Correct CL reads
		if self.clGenome != "" and os.path.exists(self.clGenome) and os.path.getsize(self.clReads) >0:
			self.refine_organelle("CL", self.clGenome, self.clReads)

		# Correct reads assigned to genes
		if self.assignedReads != "" and os.path.exists(self.assignedReads) and os.path.getsize(self.assignedReads) >0:
			self.refine_Assigned()

		# Correct reads aligned to the reference genome but not assigned to genes
		if self.notassignedReads != "" and os.path.exists(self.notassignedReads) and os.path.getsize(self.notassignedReads) >0:
			self.refine_NotAssigned()

	def refine_organelle(self, organelle, orgGenomePath, orgReadPath):

		"""
		This function is used to correct the polyA tails of CCS reads
		which are assigned to organelles (mitchondrial and choloroplast)
		"""

		# Set output directory
		outDir = ""
		if organelle == "MT":
			print("\nStart to correct polyA tails of mitchondrial reads...")

			# Create output directory
			mtDir = os.path.abspath(os.path.join(self.refinePolyaDir, 'MT'))
			Misc.Misc.checkDir(mtDir, 'MT') # check and create output dir
			outDir = mtDir

		elif organelle == "CL":
			print("\nStart to correct polyA tails of chloroplast reads...")

			# Create output directory
			clDir = os.path.abspath(os.path.join(self.refinePolyaDir, 'CL'))
			Misc.Misc.checkDir(clDir, 'CL') # check and create output dir
			outDir = clDir

		# Build index of organelle genomes
		indexPrefix = os.path.abspath(os.path.join(self.refinePolyaDir, outDir, os.path.basename(orgGenomePath)))
		indexLog = os.path.join(outDir, 'index.log')

		cmd = [self.bwaPath,'index -p', indexPrefix, orgGenomePath, '&>' + indexLog]
		cmd = ' '.join([str(x) for x in cmd])
		if organelle == "MT":
			print("\nIndex mitchondiral genome:")
		elif organelle == "CL":
			print("\nIndex chloroplast genome:")
		print(cmd)
		os.system(cmd)

		# Get draft polyA tail
		polyAPath = os.path.abspath(os.path.join( outDir, 'polyA.txt'))
		cmd = ['fisher.pl', orgReadPath, self.polyaPath, "| awk '{print $1\"\\t\"$NF}' >", polyAPath]
		cmd = ' '.join([str(x) for x in cmd])
		print('\nExtract draft polyA:')
		print(cmd)
		os.system(cmd)

		# Extract mt reads (35 nt + 10 nt)
		orgReadFaPath = os.path.abspath(os.path.join(outDir, 'refine_reads.fasta'))
		cmd = ['fisher.pl', orgReadPath, self.polyaPath, "| awk '{print \">\"$1\"\\n\"$2}' >", orgReadFaPath]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nExtract mt reads (35 nt + 10 nt):")
		print(cmd)
		os.system(cmd)

		# Map mt reads (35 nt + 10 nt) to mt genomes
		out_bam = os.path.abspath(os.path.join(outDir,'refine_reads.bam'))
		err_txt = os.path.abspath(os.path.join(outDir,'align.log'))
		out = '| ' + self.samtoolsPath + ' view -bS -F3844 > ' + out_bam
		err = '2> ' + err_txt
		cmd = [self.bwaPath,'mem', indexPrefix, '-T 20 -k11 -W10 -r10 -A1 -B1 -O1 -E1 -L0', orgReadFaPath, err, out]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nAlign:")
		print(cmd)
		os.system(cmd)

		# Extract 3' soft clip sequence
		print("\nExtract polyA:")
		self.extractPolyA(out_bam,polyAPath,outDir)

	def refine_Assigned(self):

		# Dealing with assigned reads, Create dir
		print("\nStart to correct polyA of reads assigned to genes...")
		refineReadDir = os.path.abspath(os.path.join(self.assignedDir,'refineReads'))
		refineAlignDir = os.path.abspath(os.path.join(self.assignedDir,'refineAlign'))
		Misc.Misc.checkDir(refineReadDir, 'refineReads')
		Misc.Misc.checkDir(refineAlignDir, 'refineAlign')
		outDir = self.assignedDir

		# Get draft polyA tail
		polyAPath = os.path.abspath(os.path.join( self.assignedDir, 'polyA.txt'))
		cmd = ['fisher.pl', self.assignedReads , self.polyaPath, "| awk '{print $1\"\\t\"$NF}' >", polyAPath]
		cmd = ' '.join([str(x) for x in cmd])
		print('\nExtract draft polyA:')
		print(cmd)
		os.system(cmd)

		# Extract assigned reads (35 nt + 10 nt)
		# Get assigned gene of each read
		assigned_dict = {}
		genes = {}
		with open(self.assignedReads,'r') as AR:
			for line in AR:
				(read_name, gene) = line.strip('\n').split('\t')
				assigned_dict[read_name] = gene
				genes[gene] = 1

		# Prepare reads
		refineReads_dict = {}
		with open(self.polyaPath, 'r') as PA:
			for line in PA:
				(read_name, seq10, seq30, seqpa, pa) = line.strip('\n').split('\t')

				if read_name in assigned_dict:
					gene = assigned_dict[read_name]
				else:
					continue

				if gene not in refineReads_dict:
					refineReads_dict[ gene ] = []
				refineReads_dict[ gene].append( [read_name, seq10])

		for gene, reads in refineReads_dict.items():
			ofilePath = os.path.abspath(os.path.join(refineReadDir, gene))
			ofile = open(ofilePath, 'w')
			for read in reads:
				read_name = read[0]
				read_seq10 = read[1]
				record = ['>' + read_name, read_seq10]
				record = '\n'.join([str(x) for x in record])
				ofile.write(record + '\n')
			ofile.close()

		# Map assigned reads (35 nt + 10 nt) to genes
		refineAlignBam = os.path.abspath(os.path.join(self.assignedDir,'refineAlign.bam'))
		refineAlignCmdFile = os.path.abspath(os.path.join(self.assignedDir,"refineAlignCmd.sh"))
		refineSortCmdFile = os.path.abspath(os.path.join(self.assignedDir,"refineSortCmd.sh"))
		alignCmd = open(refineAlignCmdFile, 'w')
		sortCmd = open(refineSortCmdFile, 'w')
		out_bam_list = []
		for gene in genes.keys():
			out_bam = os.path.abspath(os.path.join(refineAlignDir, gene + '.bam'))
			err_txt = os.path.abspath(os.path.join(refineAlignDir, gene + '.align.log'))
			out = '| ' + self.samtoolsPath + ' view -bS -F3844 > ' + out_bam
			err = '2> ' + err_txt
			indexPrefix = os.path.abspath(os.path.join(self.geneIndex, gene))
			faIn = os.path.abspath(os.path.join(refineReadDir, gene))
			# align cmd
			cmd = [self.bwaPath,'mem', indexPrefix, '-T 20 -k11 -W10 -r10 -A1 -B1 -O1 -E1 -L0', faIn, err, out]
			cmd = ' '.join([str(x) for x in cmd])
			alignCmd.write(cmd + '\n')

			# sort cmd
			out_srt_bam = os.path.abspath(os.path.join(refineAlignDir, gene + '.srt.bam'))
			cmd = [self.samtoolsPath ,'sort', '-o' , out_srt_bam, out_bam]
			cmd = ' '.join([str(x) for x in cmd])
			out_bam_list.append(out_srt_bam)
			sortCmd.write(cmd + '\n')

		alignCmd.close()
		sortCmd.close()

		# Align and sort
		cmd = ['ParaFly', '-c', refineAlignCmdFile, '-CPU', self.thread]
		cmd = ' '.join([str(x) for x in cmd])
		os.system(cmd)

		cmd = ['ParaFly', '-c', refineSortCmdFile, '-CPU', self.thread]
		cmd = ' '.join([str(x) for x in cmd])
		os.system(cmd)

		# merge bam files
		bamChunks = self.chunks(out_bam_list, 500)
		bamChunksMergedBam = []
		i = 0
		for bamChunk in bamChunks:
			i = i + 1
			tmp = os.path.abspath(os.path.join(self.assignedDir, refineAlignBam + str(i)))
			bamChunksMergedBam.append(tmp)

			cmd = [self.samtoolsPath, 'merge', '-f', tmp]
			cmd.extend(bamChunk)
			cmd = ' '.join([str(x) for x in cmd])
			os.system(cmd)
		cmd = [ self.samtoolsPath, 'merge', '-f', refineAlignBam ]
		cmd.extend(bamChunksMergedBam)
		cmd = ' '.join([str(x) for x in cmd])
		os.system(cmd)

		# clean
		cmd = ['rm']
		cmd.extend(bamChunksMergedBam)
		cmd = ' '.join([str(x) for x in cmd])
		os.system(cmd)

		# Extract 3' soft clip sequence
		print("\nExtract polyA:")
		self.extractPolyA(refineAlignBam,polyAPath,outDir)


	def refine_NotAssigned(self):

		"""
		This function is used to corret the polyA tails of reads
		which are aligned to the reference genome, but not assigned to any genes.
		"""

		print('\nStart to correct polyA tail of not-assigned reads...')

		# Create dir
		refineReadDir = os.path.abspath(os.path.join(self.notAssignedDir,'refineReads'))
		refineAlignDir = os.path.abspath(os.path.join(self.notAssignedDir,'refineAlign'))
		Misc.Misc.checkDir(refineReadDir, 'refineReads')
		Misc.Misc.checkDir(refineAlignDir, 'refineAlign')
		outDir = self.notAssignedDir

		# Cluster Reads
		# get bam header
		bamHeader = os.path.abspath(os.path.join(self.notAssignedDir, 'bam_header.txt'))
		cmd = [self.samtoolsPath, 'view', '-H', self.bamFile, '>', bamHeader]
		cmd = ' '.join([str(x) for x in cmd])
		print('\nGet bam header:')
		print(cmd)
		os.system(cmd)

		# Get draft polyA tail
		polyAPath = os.path.abspath(os.path.join( self.notAssignedDir, 'polyA.txt'))
		cmd = ['fisher.pl', self.notassignedReads , self.polyaPath, "| awk '{print $1\"\\t\"$NF}' >", polyAPath]
		cmd = ' '.join([str(x) for x in cmd])
		print('\nExtract draft polyA:')
		print(cmd)
		os.system(cmd)

		# not assigned reads (BAM)
		notAssignedBamPath = os.path.abspath(os.path.join(self.notAssignedDir, 'NotAssigned.bam'))
		cmd = [ self.samtoolsPath, 'view', self.bamFile, '| fisher.pl', self.notassignedReads, '- | cat', bamHeader,'- | ', self.samtoolsPath, 'view', '-bS', '>' + notAssignedBamPath ]
		cmd = ' '.join( [str(x) for x in cmd] )
		print(cmd)
		os.system(cmd)

		# not assigned reads (BED)
		notAssignedBedPath = os.path.abspath(os.path.join(self.notAssignedDir, 'NotAssigned.bed'))
		cmd = [self.bedtoolsPath, 'bamtobed', '-i', notAssignedBamPath, '>' + notAssignedBedPath]
		cmd = ' '.join([ str(x) for x in cmd])
		print("\nConvert BAM to BED:")
		print(cmd)
		os.system(cmd)

		# calculate intervals
		notAssignedBedFilterPath=os.path.abspath(os.path.join(self.notAssignedDir,'NotAssigned.filtered.bed'))
		cmd = [ 'awk ', "'$3 - $2 < 100000'", notAssignedBedPath, '>' + notAssignedBedFilterPath ]
		cmd = ' '.join([ str(x) for x in cmd])
		print('\nFilter reads:')
		print(cmd)
		os.system(cmd)

		# cluster
		clusterPath = os.path.abspath( os.path.join(self.notAssignedDir, 'NotAssigned.cluster') )
		cmd=[self.bedtoolsPath, 'cluster', '-s', '-i', notAssignedBedFilterPath, '-d 200', '>' + clusterPath]
		cmd=' '.join([ str(x) for x in cmd])
		print('\nCluster reads:')
		print(cmd)
		os.system(cmd)

		# filter cluster (CCS >= 2)
		filteredClusterPath=os.path.abspath(os.path.join(self.notAssignedDir,'NotAssigned.cluster.filtered') )
		cmd = [ 'cat', clusterPath, "| awk '{print $NF}' | uniq -c | awk '$1 >= 2 {print $2}' | fisher.pl -fc 7 -", clusterPath, '>' + filteredClusterPath]
		cmd = ' '.join([ str(x) for x in cmd])
		print('\nCluster reads (CCS >= 2):')
		print(cmd)
		os.system(cmd)

		# split reads by cluster
		clusterDir = os.path.abspath(os.path.join(self.notAssignedDir, 'clusters'))
		Misc.Misc.checkDir(clusterDir, 'clusters')
		record = ''
		last_cluster_name = ''
		clusters = {}
		cluster_dict = {}
		with open(filteredClusterPath, 'r') as f:
			for line in f:
				(chrom, start, end, read_name, score, strand, cluster_name) = line.strip('\n').split('\t')
				clusters[cluster_name] = 1
				cluster_dict[read_name] = cluster_name
				if last_cluster_name == '':
					last_cluster_name = cluster_name
				if cluster_name != last_cluster_name:
					ofilePath = os.path.abspath(os.path.join(clusterDir, last_cluster_name + '.bed'))
					of = open(ofilePath, 'w')
					of.write(record)
					of.close()
					record = line
				else:
					record = record + line
				last_cluster_name = cluster_name
		ofilePath = os.path.abspath(os.path.join(clusterDir, last_cluster_name + '.bed'))
		of = open(ofilePath, 'w')
		of.write(record)
		of.close()

		# calculate interval of clusters
		print('\nCalculate interval of clusters')
		mergeClusterDir = os.path.abspath(os.path.join(self.notAssignedDir, 'merge'))
		Misc.Misc.checkDir(mergeClusterDir, 'merge')
		intervals = []
		for cluster in clusters:
			infile = os.path.abspath(os.path.join(clusterDir, cluster + '.bed'))
			intervalPath = os.path.abspath(os.path.join(mergeClusterDir, cluster + '.bed'))
			intervals.append(intervalPath)
			cmd =[ self.bedtoolsPath , 'merge', '-d', 200, '-c', '6,7', '-o', 'distinct,distinct', '-s', '-i', infile, '>' + intervalPath ]
			cmd = ' '.join([ str(x) for x in cmd])
			os.system(cmd)

		# combine cluster intervals
		print("\nCollect intervals")
		allIntervalPath = os.path.abspath(os.path.join(self.notAssignedDir, 'cluster_interval.bed'))
		allIntervalFile = open(allIntervalPath, 'w')
		for interval in intervals:
			with open(interval, 'r') as IV:
				for line in IV:
					allIntervalFile.write(line)
			IV.close()
		allIntervalFile.close()

		# Calculate genome size
		chromSizePath = os.path.abspath(os.path.join(self.notAssignedDir,"chrom_size.txt"))
		chromSizeFile = open(chromSizePath, 'w')
		with pysam.FastxFile(self.genomeFasta) as fa_in:
			for cnt, read in enumerate(islice(fa_in,None)):
				chrom_size = len(read.sequence)
				chromSizeFile.write( read.name + '\t' + str(chrom_size) + '\n')
		chromSizeFile.close()

		# Extend bed
		allIntervalExtendPath=os.path.abspath(os.path.join(self.notAssignedDir,'cluster_interval.extend.bed'))
		cmd = [self.bedtoolsPath, 'slop', '-l 2000', '-r 2000', '-s', '-i', allIntervalPath, '-g', chromSizePath, '>' + allIntervalExtendPath]
		cmd = ' '.join( [str(x) for x in cmd] )
		os.system(cmd)

		tmp = allIntervalExtendPath + '.tmp'
		cmd = ["awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$5,0,$4}'", allIntervalExtendPath ,'>' + tmp, '&& mv', tmp , allIntervalExtendPath]
		cmd = ' '.join( [str(x) for x in cmd] )
		print(cmd)
		os.system(cmd)

		# Get dna sequence of extended intervals
		intervalFaPath = os.path.abspath(os.path.join(self.notAssignedDir,'cluster_interval.extend.fa'))
		cmd = [self.bedtoolsPath, 'getfasta' , '-name', '-fi', self.genomeFasta, '-bed', allIntervalExtendPath, '>' + intervalFaPath]
		cmd = ' '.join( [str(x) for x in cmd] )
		os.system(cmd)

		# Build index for each cluster
		geneDir = os.path.abspath(os.path.join(self.notAssignedDir, 'gene'))
		indexDir = os.path.abspath(os.path.join(self.notAssignedDir, 'index'))
		Misc.Misc.checkDir(geneDir, 'gene')
		Misc.Misc.checkDir(indexDir, 'index')

		indexCmdPath = os.path.abspath(os.path.join(self.notAssignedDir, 'indexCmd.sh'))
		indexCmdFile = open(indexCmdPath, 'w')
		with pysam.FastxFile(intervalFaPath) as fa_file:
			for cnt, e in enumerate(islice(fa_file,None)):
				ofPath = os.path.abspath(os.path.join(geneDir, e.name))
				of = open(ofPath, 'w')
				of.write(str(e))
				of.close()

				idxPath = os.path.abspath(os.path.join(indexDir, e.name))
				indexLog = os.path.abspath(os.path.join(indexDir, e.name +'_index.log'))
				cmd = [self.bwaPath,'index -p', idxPath, ofPath, '&>', indexLog]
				cmd = ' '.join( [str(x) for x in cmd] )
				indexCmdFile.write(cmd + '\n')
		indexCmdFile.close()

		# Run build index
		cmd = ['ParaFly', '-c', indexCmdPath, '-CPU', self.thread]
		cmd = ' '.join([str(x) for x in cmd])
		print('\nIndex cluster intervals: ' + cmd)
		os.system(cmd)

		# Prepare reads
		refineReads_dict = {}
		with open(self.polyaPath, 'r') as PA:
			for line in PA:
				(read_name, seq10, seq30, seqpa, pa) = line.strip('\n').split('\t')
				gene = ''
				if read_name in cluster_dict:
					gene = cluster_dict[read_name]
				else:
					continue

				if gene not in refineReads_dict:
					refineReads_dict[ gene ] = []
				refineReads_dict[ gene].append( [read_name, seq10])

		genes = {}
		for gene, reads in refineReads_dict.items():
			genes[gene] = 1
			ofilePath = os.path.abspath(os.path.join(refineReadDir, gene))
			ofile = open(ofilePath, 'w')
			for read in reads:
				read_name = read[0]
				read_seq10 = read[1]
				record = ['>' + read_name, read_seq10]
				record = '\n'.join([str(x) for x in record])
				ofile.write(record + '\n')
			ofile.close()

		# Map notassigned reads (35 nt + 10 nt) to clusters
		refineAlignBam = os.path.abspath(os.path.join(self.notAssignedDir,'refineAlign.bam'))
		refineAlignCmdFile = os.path.abspath(os.path.join(self.notAssignedDir,"refineAlignCmd.sh"))
		refineSortCmdFile = os.path.abspath(os.path.join(self.notAssignedDir,"refineSortCmd.sh"))
		alignCmd = open(refineAlignCmdFile, 'w')
		sortCmd = open(refineSortCmdFile, 'w')
		out_bam_list = []
		for gene in genes.keys():
			out_bam = os.path.abspath(os.path.join(refineAlignDir, gene + '.bam'))
			err_txt = os.path.abspath(os.path.join(refineAlignDir, gene + '.align.log'))
			out = '| ' + self.samtoolsPath + ' view -bS -F3844 > ' + out_bam
			err = '2> ' + err_txt
			indexPrefix = os.path.abspath(os.path.join(indexDir, gene))
			faIn = os.path.abspath(os.path.join(refineReadDir, gene))

			# align cmd
			cmd = [self.bwaPath,'mem', indexPrefix, '-T 20 -k11 -W10 -r10 -A1 -B1 -O1 -E1 -L0', faIn, err, out]
			cmd = ' '.join([str(x) for x in cmd])
			alignCmd.write(cmd + '\n')

			# sort cmd
			out_srt_bam = os.path.abspath(os.path.join(refineAlignDir, gene + '.srt.bam'))
			cmd = [self.samtoolsPath ,'sort', '-o' , out_srt_bam, out_bam]
			cmd = ' '.join([str(x) for x in cmd])
			out_bam_list.append(out_srt_bam)
			sortCmd.write(cmd + '\n')

		alignCmd.close()
		sortCmd.close()

		# Align and sort
		cmd = ['ParaFly', '-c', refineAlignCmdFile, '-CPU', self.thread]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nNot Assigned, Align: " + cmd)
		os.system(cmd)

		cmd = ['ParaFly', '-c', refineSortCmdFile, '-CPU', self.thread]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nNot Assigned, Sort: " + cmd)
		os.system(cmd)

		# Merge bam files
		bamChunks = self.chunks(out_bam_list, 500)
		bamChunksMergedBam = []
		i = 0
		for bamChunk in bamChunks:
			i = i + 1
			tmp = os.path.abspath(os.path.join(self.notAssignedDir, refineAlignBam + str(i)))
			bamChunksMergedBam.append(tmp)

			cmd = [self.samtoolsPath, 'merge', '-f', tmp]
			cmd.extend(bamChunk)
			cmd = ' '.join([str(x) for x in cmd])
			os.system(cmd)
		cmd = [ self.samtoolsPath, 'merge', '-f', refineAlignBam ]
		cmd.extend(bamChunksMergedBam)
		cmd = ' '.join([str(x) for x in cmd])
		os.system(cmd)

		# clean
		cmd = ['rm']
		cmd.extend(bamChunksMergedBam)
		cmd = ' '.join([str(x) for x in cmd])
		os.system(cmd)

		# Extract 3' soft clip sequence
		print("\nExtract polyA...")
		self.extractPolyA(refineAlignBam,polyAPath,outDir)

	def extractPolyA(self,bamfile,draftfile,outDir):

		outfile = os.path.abspath(os.path.join( outDir, os.path.basename(outDir)+ '.refine.polyA.txt'))
		fout = open(outfile, 'w')

		dic_pa = {}
		with open(draftfile, 'r') as FP:
			for line in FP:
				(ID,paSeq) = line.strip().split("\t")
				dic_pa[ID] = paSeq

		prefix = bamfile.rstrip('.bam')
		pairwiseFile =  prefix + '.sam2pairwise.out'
		print("\nConvert bam format to pairwise format:")
		cmd = [self.samtoolsPath,'view',bamfile, '|', self.sam2pairwisePath, '>',pairwiseFile]
		cmd = ' '.join([ str(x) for x in cmd])
		print(cmd)
		os.system(cmd)

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
							read = read.replace('-','')
							read = read.replace('.','')
							break
						else:
							continue
				else:
					for i in ref:
						if i == "N":
							n += 1
						else:
							break
					for j in range(n,len(ref),1):
						if j == len(ref) - 1:
							read = read.replace('-','')
							read = read.replace('.','')
							read = self.revComp(read)
							tail = read
						elif read[j] != "T" and ref[j] != "T" and read[j] == ref[j]:
							tail = read[:j]
							tail = tail.replace('-','')
							tail = tail.replace('.','')
							tail = self.revComp(tail)
							read = read.replace('-','')
							read = read.replace('.','')
							read = self.revComp(read)
							break
						else:
							continue
			
			#while "G" in three_soft_clip_seq[:2] or "C" in three_soft_clip_seq[:2]:
			#	three_soft_clip_seq = three_soft_clip_seq[1:]
			#	if three_soft_clip_seq == "":
			#		break
				polya_seq = tail + dic_pa[ID][10:]
				if polya_seq != "":
					fout.write(ID + "\t" + polya_seq + "\n")
		fout.close()

	def chunks(self,lst,n):
		"""Yield successive n-sized chunks from lst."""
		for i in range(0, len(lst), n):
			yield lst[i:i + n]

	def revComp(self,s):
		revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
		return ''.join([revCDict[x] for x in s])[::-1]
