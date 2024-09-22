#!/usr/bin/env python
import argparse
import PAIsoSeqMaster

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Process PAIso-Seq Sequencing Data')
	subparser = parser.add_subparsers(dest="command")

	# demultiplex and clean Raw CCS reads
	parser_demultiplex = subparser.add_parser('preprocess', help="Demultiplex and clean CCS reads")
	parser_demultiplex.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	# align/map Clean CCS reads to reference genomes
	parser_aligntogenome = subparser.add_parser('aligntoGenome', help="Align CCS reads to reference genomes")
	parser_aligntogenome.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	# assign CCS reads to genes
	parser_assign = subparser.add_parser('assign', help="Assign CCS reads to genes")
	parser_assign.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	# align bad polyA reads to cDNA
	parser_aligntotranscript = subparser.add_parser('aligntoTranscript', help="Align bad polyA ccs reads to genes")
	parser_aligntotranscript.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	# classification of each CCS reads
	parser_classify = subparser.add_parser('classify', help="Classification of each CCS reads")
	parser_classify.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)
	# define draft poly(A) tail
	parser_draft = subparser.add_parser('draft', help="Define draft poly(A) tail")
	parser_draft.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	# refine poly(A) tail
	parser_refine = subparser.add_parser('refine',help="Refine poly(A) tail")
	parser_refine.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	# collect result
	parser_result = subparser.add_parser('result',help="Collect results")
	parser_result.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	# All
	parser_all = subparser.add_parser('all', help="Run Complete PAIso-Seq Pipeline")
	parser_all.add_argument('-p', '--parameters', help="Path To Experiment parameter.yaml", required=True)

	params = parser.parse_args()
	command = params.command
	parameters = params.parameters

	if command == 'all':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.preprocess()
		master.aligntoGenome()
		master.assign()
		master.aligntoTranscript()
		master.classify()
		master.draft()
		master.refine()
		master.result()
	elif command == 'preprocess':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.preprocess()
	elif command == 'aligntoGenome':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.aligntoGenome()
	elif command == 'assign':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.assign()
	elif command == 'aligntoTranscript':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.aligntoTranscript()
	elif command == 'classify':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.classify()
	elif command == 'draft':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.draft()
	elif command == 'refine':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.refine()
	elif command == 'result':
		master = PAIsoSeqMaster.PAIsoSeqMaster(parameterYamlFile=parameters)
		master.results()
	else:
		print('Please Specify Correct Command!')
