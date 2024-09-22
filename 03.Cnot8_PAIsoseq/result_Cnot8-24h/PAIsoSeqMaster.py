#!/usr/bin/env python
import os
import sys
from Source import Misc
from Source import Parameters
from Source import Preprocess
from Source import AligntoGenome
from Source import Assign
from Source import AligntoTranscript
from Source import Classify
from Source import Draft
from Source import Refine
from Source import Results

class PAIsoSeqMaster():

	def __init__(self, parameterYamlFile):
		print("------------------")
		print("PAIso-Seq Pipeline")
		print("------------------")

		# Generate Parameters Object
		self.parameters = Parameters.Parameters(parameterYamlFile)

		# Check outputDir
		expOutDir = self.parameters.getParametersBysubKey('experiment', 'outputDir')
		Misc.Misc.checkDir(expOutDir, 'OutputDir')


	def getParameters(self):
		return self.parameters.getParameters()

	def preprocess(self):
		preprocess = Preprocess.Preprocess(self.parameters)
		preprocess.process()

	def aligntoGenome(self):
		aligntoGenome = AligntoGenome.AligntoGenome(self.parameters)
		aligntoGenome.aligntoGenome()

	def assign(self):
		assign = Assign.Assign(self.parameters)
		assign.assign()

	def aligntoTranscript(self):
		aligntoTranscript = AligntoTranscript.AligntoTranscript(self.parameters)
		aligntoTranscript.aligntoTranscript()

	def classify(self):
		classify = Classify.Classify(self.parameters)
		classify.classify()

	def draft(self):
		draft = Draft.Draft(self.parameters)
		draft.draft()

	def refine(self):
		refine = Refine.Refine(self.parameters)
		refine.refine()

	def results(self):
		results = Results.Results(self.parameters)
		results.results()
