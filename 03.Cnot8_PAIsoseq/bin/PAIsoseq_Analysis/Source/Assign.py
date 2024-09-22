#!/usr/bin/env python
import os
import subprocess
from Source import Misc

class Assign():

	def __init__(self, parameters):

		print("\n---------------------")
		print("Assign reads to genes")
		print("---------------------")

		# Output directory
		outputDir = parameters.getParametersBysubKey('experiment', 'outputDir')

		# Input bam file
		alignDir = os.path.join(outputDir, 'AligntoGenomeDir')
		#self.bamPath = os.path.join(alignDir,'clean_CCS.sorted.bam')
		self.bamPath = os.path.join(alignDir,'clean_CCS.sorted.all_mapped.bam')
		Misc.Misc.checkFile(self.bamPath, 'Alignment file')

        # Generate new dir for writing output files
		self.assignDir = os.path.join(outputDir, 'AssignDir')
		Misc.Misc.checkDir(self.assignDir, 'AssignDir')

        # Check if software is accressible
		self.featureCountPath = parameters.getParametersBysubKey('software', 'featureCounts')
		Misc.Misc.checkProg(self.featureCountPath, 'featureCounts')

        # Check GTF File
		self.gtfFile = parameters.getParametersBysubKey('experiment', 'annotationGTF')
		Misc.Misc.checkFile(self.gtfFile, 'Annotation GTF')

	def assign(self):
		self.quantifyFeatCounts()

	def quantifyFeatCounts(self):
		featCountPath = os.path.join(self.assignDir,'clean_CCS.fcounts')
		cmd = [self.featureCountPath, '-L', '-g', 'gene_id', '-t', 'exon', '-s', str(1), '-R', 'CORE', '-a', self.gtfFile, '-o', featCountPath, self.bamPath]
		cmd = ' '.join([str(x) for x in cmd])
		print(cmd)
		os.system(cmd)
