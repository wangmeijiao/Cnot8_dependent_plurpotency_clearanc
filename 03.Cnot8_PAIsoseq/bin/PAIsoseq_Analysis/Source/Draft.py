#!/usr/bin/env python
import os
import subprocess
from Source import Misc

class Draft():

	def __init__(self, parameters):

		print("\n-------------------------")
		print("Define draft poly(A) tail")
		print("-------------------------")

		# Output directory
		outputDir = parameters.getParametersBysubKey('experiment', 'outputDir')

		# Input file
		aligntoTranscriptDir = os.path.join(outputDir, 'AligntoTranscriptDir')
		self.badPolyaPath = os.path.join(aligntoTranscriptDir, 'clean_CCS.sorted.all_mapped.badPolyA.refine.txt')
		self.badPolyaPath = os.path.abspath(self.badPolyaPath)
		Misc.Misc.checkFile(self.badPolyaPath, '3-soft-clip containing nonA bases')

		# Generate new dir for writing output files
		self.draftPolyaDir= os.path.join(outputDir, 'DraftPolyaDir')
		Misc.Misc.checkDir(self.draftPolyaDir, 'DraftPolyaDir')

	def draft(self):
		# link data
		badPolyaPath = os.path.join(self.draftPolyaDir, os.path.basename(self.badPolyaPath) )
		if os.path.exists(badPolyaPath) == False:
			cmd = ['ln -s', self.badPolyaPath, self.draftPolyaDir]
			cmd = ' '.join([str(x) for x in cmd])
			os.system(cmd)
		self.badPolyaPath = os.path.abspath(badPolyaPath)

		# Define draft poly(A) tail
		out_txt = os.path.join(self.draftPolyaDir, 'clean_CCS_sorted_all_mapped_badPolyA_refine_draft.txt')
		err_txt = os.path.join(self.draftPolyaDir, 'clean_CCS_sorted_all_mapped_badPolyA_refine_draft.err')
		out = '1>' + out_txt
		err = '2> ' + err_txt
		cmd = ['getPolyA.py', self.badPolyaPath, 10, out, err]
		cmd = ' '.join( [ str(x) for x in cmd])
		print(cmd)
		os.system(cmd)
