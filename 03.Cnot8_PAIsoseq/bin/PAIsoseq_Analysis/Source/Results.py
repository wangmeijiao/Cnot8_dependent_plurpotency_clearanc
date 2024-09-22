import os
#import pandas as pd
from Source import Misc

class Results():

	def __init__(self, parameters):
		print('\n---------------')
		print('Collect Results')
		print('---------------')

		outputDir = parameters.getParametersBysubKey('experiment', 'outputDir')

		# Create output directory
		self.resultDir = os.path.abspath(os.path.join(outputDir, 'ResultDir'))
		Misc.Misc.checkDir(self.resultDir, 'ResultDir')

		# Results
		self.nicepolyaResult = os.path.abspath(os.path.join(outputDir, 'AligntoGenomeDir', 'clean_CCS.sorted.all_mapped.nicePolyA.txt'))
		self.badpolyaResult = os.path.abspath(os.path.join(outputDir, 'AligntoTranscriptDir', 'badPolya_3sc_window0.txt'))
		self.mtResult = os.path.abspath(os.path.join(outputDir, 'RefinePolyaDir', 'MT', 'MT.refine.polyA.txt'))
		self.clResult = os.path.abspath(os.path.join(outputDir, 'RefinePolyaDir', 'CL', 'CL.refine.polyA.txt'))
		self.assignedResult = os.path.abspath(os.path.join(outputDir, 'RefinePolyaDir', 'Assigned', 'Assigned.refine.polyA.txt'))
		self.notAssignedResult = os.path.abspath(os.path.join(outputDir, 'RefinePolyaDir', 'NotAssigned', 'NotAssigned.refine.polyA.txt'))
		self.annotationPath = os.path.abspath(os.path.join(outputDir, 'AssignDir', 'clean_CCS.sorted.all_mapped.bam.featureCounts'))

	def results(self):
		polyAPath = os.path.abspath(os.path.join(self.resultDir,'all.polyA_seq.txt'))
		cmd = ['cat']

		if os.path.exists(self.mtResult):
			cmd.extend( [self.mtResult] )
			print("\nPolyA of mitchondrial reads: " + self.mtResult)
		else:
			print(self.mtResult + "not exist!")

		if os.path.exists(self.clResult):
			cmd.extend( [self.clResult] )
			print("\nPolyA of chloroplast reads: " + self.clResult)
		else:
			print(self.clResult + "not exist!")

		if os.path.exists(self.assignedResult):
			cmd.extend( [self.assignedResult] )
			print("\nPolyA of assigned reads: " + self.assignedResult)
		else:
			print(self.assignedResult + "not exist!")

		if os.path.exists(self.notAssignedResult):
			cmd.extend( [self.notAssignedResult] )
			print("\nPolyA of not-assigned reads: " + self.notAssignedResult)
		else:
			print(self.notAssignedResult + "not exist!")

		cmd.extend( ["| awk 'NF == 2 {print $1\"\\t\"$2}' |"])
		cmd.extend( ["cat - "])

		if os.path.exists(self.badpolyaResult):
			cmd.extend( [self.badpolyaResult])
			print("\nFiltered bad polyA reads: " + self.badpolyaResult)

		if os.path.exists(self.nicepolyaResult):
			cmd.extend( [self.nicepolyaResult] )
			print("\nNice polyA reads: " + self.nicepolyaResult)

		cmd.extend( [ '>' + polyAPath] )
		cmd = ' '.join([str(x) for x in cmd])
		print("\nCombine results:")
		print(cmd)
		os.system(cmd)

		polyATailPath = os.path.abspath(os.path.join(self.resultDir,'all.polyA_tail_info.txt'))
		cmd = ['formatPolyAResult.py', polyAPath, self.annotationPath, '>' + polyATailPath]
		cmd = ' '.join([str(x) for x in cmd])
		print("\nFormat results:")
		print(cmd)
		os.system(cmd)
		print("\nResult: " + polyATailPath)
