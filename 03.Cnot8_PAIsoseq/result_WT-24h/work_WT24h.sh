

#step 0 initiate environment
conda activate pacbio
export PATH=$PATH:/home/mjwang/pwd/Quan-Cnot8-polyA/02.result_WT-24h/Script/

#Step 1: Demultiplex and clean of CCS reads # preprocess
./PAIsoSeqAnalysis.py preprocess -p parameters_138.yaml > log.step01 2>&1 &


#Step 2: Align CCS reads to reference genomes 
./PAIsoSeqAnalysis.py aligntoGenome -p parameters_138.yaml > log.step02 2>&1 &


#Step 3: Assign CCS reads to genes and clusters
./PAIsoSeqAnalysis.py assign -p parameters_138.yaml > log.step03 2>&1 &

#Step 4: Align to transcript (cDNA)
./PAIsoSeqAnalysis.py aligntoTranscript -p parameters_138.yaml > log.step04 2>&1 &


#Step 5: Classification of each CCS reads
./PAIsoSeqAnalysis.py classify -p parameters_138.yaml > log.step05 2>&1 &


#Step 6: Define a draft poly(A) tail
./PAIsoSeqAnalysis.py draft -p parameters_138.yaml > log.step06 2>&1 &

#Step 7: Refine polyA(A) tail
./PAIsoSeqAnalysis.py refine -p parameters_138.yaml > log.step07 2>&1 &


#Step 8: Collect result
./PAIsoSeqAnalysis.py result -p parameters_138.yaml > log.step08 2>&1 &




########try in one step#####
conda activate pacbio
export PATH=$PATH:/home/mjwang/pwd/Quan-Cnot8-polyA/02.result_WT-24h/Script/
./PAIsoSeqAnalysis.py all -p parameters_138.yaml > log.runall.new 2>&1 &
#./PAIsoSeqAnalysis.py all -p parameters_138.yaml > log.runall 2>&1 &

#add the step08 because the all pipe failure at this step, corrected at Oct.8
#./PAIsoSeqAnalysis.py result -p parameters_138.yaml > log.step08 2>&1 &



