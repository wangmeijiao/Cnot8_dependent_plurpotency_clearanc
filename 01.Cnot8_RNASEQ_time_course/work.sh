

cat samples.all.txt |while read a b; do echo $a $b;ln -s ../03.mapping_tophat2/${b}_tophat2/accepted_hits.bam $b.accepted_hits.bam ; done

cat samples.all.txt |while read a b; do echo -ne $b".accepted_hits.bam "; done

featureCounts -T 12 -B -C -O -p -a mm9.refSeq.txt.gtf -o out.ftcount.txt Cnot8-2-WT-ESC_rep1.accepted_hits.bam Cnot8-4-WT-ESC_rep2.accepted_hits.bam Cnot8-5-KO-ESC_rep1.accepted_hits.bam Cnot8-8-KO-ESC_rep2.accepted_hits.bam Cnot8-2-WT-6h_rep1.accepted_hits.bam Cnot8-4-WT-6h_rep2.accepted_hits.bam Cnot8-5-KO-6h_rep1.accepted_hits.bam Cnot8-8-KO-6h_rep2.accepted_hits.bam Cnot8-2-WT-12h_rep1.accepted_hits.bam Cnot8-4-WT-12h_rep2.accepted_hits.bam Cnot8-5-KO-12h_rep1.accepted_hits.bam Cnot8-8-KO-12h_rep2.accepted_hits.bam Cnot8-2-WT-24h_rep1.accepted_hits.bam Cnot8-4-WT-24h_rep2.accepted_hits.bam Cnot8-5-KO-24h_rep1.accepted_hits.bam Cnot8-8-KO-24h_rep2.accepted_hits.bam Cnot8-2-WT-48h_rep1.accepted_hits.bam Cnot8-4-WT-48h_rep2.accepted_hits.bam Cnot8-5-KO-48h_rep1.accepted_hits.bam Cnot8-8-KO-48h_rep2.accepted_hits.bam > log.cnt 2> err.cnt &

#test method in paper
featureCounts -T 12 --primary  -a mm9.refSeq.txt.gtf -o out.ftcount.testpaper_method.txt Cnot8-2-WT-ESC_rep1.accepted_hits.bam Cnot8-4-WT-ESC_rep2.accepted_hits.bam Cnot8-5-KO-ESC_rep1.accepted_hits.bam Cnot8-8-KO-ESC_rep2.accepted_hits.bam Cnot8-2-WT-6h_rep1.accepted_hits.bam Cnot8-4-WT-6h_rep2.accepted_hits.bam Cnot8-5-KO-6h_rep1.accepted_hits.bam Cnot8-8-KO-6h_rep2.accepted_hits.bam Cnot8-2-WT-12h_rep1.accepted_hits.bam Cnot8-4-WT-12h_rep2.accepted_hits.bam Cnot8-5-KO-12h_rep1.accepted_hits.bam Cnot8-8-KO-12h_rep2.accepted_hits.bam Cnot8-2-WT-24h_rep1.accepted_hits.bam Cnot8-4-WT-24h_rep2.accepted_hits.bam Cnot8-5-KO-24h_rep1.accepted_hits.bam Cnot8-8-KO-24h_rep2.accepted_hits.bam Cnot8-2-WT-48h_rep1.accepted_hits.bam Cnot8-4-WT-48h_rep2.accepted_hits.bam Cnot8-5-KO-48h_rep1.accepted_hits.bam Cnot8-8-KO-48h_rep2.accepted_hits.bam > log.testpaper.method.cnt 2> err.test.papermethod.cnt &
##count not as pairs, correlation very good with paired


##featureCounts -T 24 -B -C -F SAF -O -p -a  mm9.fa.size.win50k.bed.SAF

cat out.ftcount.txt |perl format_ft.pl > out.ftcount.format.mat 2> err.format

grep -v Unassigned_Unmapped out.ftcount.txt.summary|grep -v Status |datamash sum 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21|sed 's/\t/\,/g' |less
#without unmapped reads



##gene all uniq no dup

awk 'NR > 2' out.ftcount.rep12.txt |awk -vOFS='\t' '{print $1,$6,$7,$8}' > out.ftcount.mat

#perl -nle 'chomp;next if($_=~/^#/) ;my @box = split/[\t ]+/;my @chrs = split/;/,$box[1];     ;my @starts=split/;/,$box[2]; my @ends = split/;/,$box[3];print "$box[0]\t$chrs[0]\t$starts[0]\t$ends[-1]\t$box[5]\t$box[6]\t$box[7]"; ' out.ftcount.rep12.txt > out.ftcount.format.mat

cat out.ftcount.rep12.txt |perl format_ft.pl > out.ftcount.format.mat 2> err.format

awk -vOFS='\t' 'NR>1{print $2,$3,$4,$1"-"$6"|"$7}' out.ftcount.format.mat > out.ftcount.format.mat.bed

#lib size
#59512913	62975656
awk -vOFS='\t' 'NR>1{print $2,$3,$4,$1,1e9*$6/($5*59512913),1e9*$7/($5*62975656)}' out.ftcount.format.mat > out.ftcount.format.mat.fpkm.bed

#check if start >= end, no
awk '$2>=$3' out.ftcount.format.mat.fpkm.bed






####
rm run.sh; cat samples.all.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${b}_cufflinks -G mm9.ucsc.knowGene.txt.gtf  ${b}.accepted_hits.bam > log.cufflinks.$b 2>&1 & \n" >> run.sh ; done

rm run_refGene.sh; cat samples.all.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${b}_cufflinks_refGene -G refGene.mm9.gtf  ${b}.accepted_hits.bam > log.cufflinks_refGene.$b 2>&1 & \n" >> run_refGene.sh ; done #genePred2gtf downloaded
rm run_refGene_reorg.sh; cat samples.all.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${b}_cufflinks_refGene_reorg -G refGene.mm9.reorganize.gtf  ${b}.accepted_hits.bam > log.cufflinks_refGene_reorg.$b 2>&1 & \n" >> run_refGene_reorg.sh ; done

rm run_refGene2.sh; cat samples.all.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${b}_cufflinks_refGene2 -G mm9.refSeq.txt.gtf  ${b}.accepted_hits.bam > log.cufflinks_refGene2.$b 2>&1 & \n" >> run_refGene2.sh ; done #ucsc table browser downloaded

rm run_refGene2.1.sh; cat samples.all.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${b}_cufflinks_refGene2.1 -G mm9.refSeq.txt.1.gtf  ${b}.accepted_hits.bam > log.cufflinks_refGene2.1.$b 2>&1 & \n" >> run_refGene2.1.sh ; done #ucsc table browser downloaded, use refSeq id as transcript id

rm run_igv.sh; cat samples.all.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${b}_cufflinks_igvGene -G mm9.igvGene.txt.gtf  ${b}.accepted_hits.bam > log.cufflinks_igvGene.$b 2>&1 & \n" >> run_igvGene.sh ; done




#rm run_protein_coding.sh; cat samples.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${b}_cufflinks_protein_coding -G Mus_musculus.GRCm38.92.rename.protein_coding.gtf  ${b}.accepted_hits.bam > log.cufflinks.protein_coding.$b 2>&1 & \n" >> run_protein_coding.sh ; done
#bash run_protein_coding.sh


cat samples.txt |while read a b; do echo $a $b;
cat ../${a}_cufflinks/genes.fpkm_tracking |perl /home/mjwang/progs/misc-tools/fpkmTools/extractFPKM.pl -log -add 1  | sort -k 4,4 > ${a}_genes.log2fpkm_tracking.sort.bed ;
done



############
cat samples.txt |while read a b; do echo $a $b;ln -s ../03.mapping_tophat2/${a}_tophat2_mcherry/accepted_hits.bam $a.accepted_hits.bam ; done

##do cufflinks
rm run.sh; cat samples.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${a}_cufflinks -G refGene.mm9Tomm10.addRFP.gtf  ${a}.accepted_hits.bam > log.cufflinks.$a 2>&1 & \n" >> run.sh ; done 

rm run_protein_coding.sh; cat samples.txt |while read a b; do echo $a $b; echo  -ne "cufflinks -p 6 -o ${a}_cufflinks_protein_coding -G Mus_musculus.GRCm38.92.rename.protein_coding.gtf  ${a}.accepted_hits.bam > log.cufflinks.protein_coding.$a 2>&1 & \n" >> run_protein_coding.sh ; done



cat samples.txt |while read a b; do echo $a $b;
cat ../${a}_cufflinks/genes.fpkm_tracking |perl /home/mjwang/progs/misc-tools/fpkmTools/extractFPKM.pl -log -add 1  | sort -k 4,4 > ${a}_genes.log2fpkm_tracking.sort.bed ;
done


##########


#2, cufflinks (-G quantify given gff only); cuffquan for cuffnorm input

cufflinks -p 20 -o SRR7695406_cuff -G mm10.ucsc.refSeq.gtf accepted_hits.bam > log.SRR7695406 2>&1 &

cat cufflink_out/isoforms.fpkm_tracking |perl /home/mjwang/progs/misc-tools/fpkmTools/extractFPKM.pl > cufflink_out/isoforms.fpkm_tracking.bed








ln -s /home/mjwang/pwd/heterochromatin_sonicate/02.analysis/domain_calling_zaret2019/zdomain_caller/mm9_zaret2019_ZDomainCaller_step1_10k_1k.10kbWindows.1kbSlide.EnrichmentScores.2-RepAverage.DomainsAfter_1.1_Cutoff.PrunedEdges.Merged.bed ./srHC_large-small.zdomain.bed

ln -s /home/mjwang/pwd/heterochromatin_sonicate/02.analysis/domain_calling_zaret2019/zdomain_caller/mm9_zaret2019_ZDomainCaller_10k_1k.euchr.bed .

#use >25% to defined a srHC-gene or euchr-gene

intersectBed -f 0.25 -a out.ftcount.format.mat.fpkm.bed -b srHC_large-small.zdomain.bed -wao > out.ftcount.format.mat.fpkm.bed.overlap_srHC

intersectBed -f 0.25 -a out.ftcount.format.mat.fpkm.bed -b mm9_zaret2019_ZDomainCaller_10k_1k.euchr.bed -wao > out.ftcount.format.mat.fpkm.bed.overlap_euchr

perl classify.pl out.ftcount.format.mat.fpkm.bed out.ftcount.format.mat.fpkm.bed.overlap_srHC out.ftcount.format.mat.fpkm.bed.overlap_euchr > result.classify 2> err.classify
grep "euch" result.classify > result.classify.euchr  #10273
grep "srHC" result.classify > result.classify.srHC  #9025

grep "few" result.classify |wcl #982
grep "both" result.classify |wcl #4027












###########
awk -vOFS='\t' '{print $4,$1,$2,$3,"+"}' mm9.fa.size.win10k.step500.bed > mm9.fa.size.win10k.step500.bed.SAF
awk -vOFS='\t' '{print $4,$1,$2,$3,"+"}' mm9.fa.size.win10k.bed > mm9.fa.size.win10k.bed.SAF
awk -vOFS='\t' '{print $4,$1,$2,$3,"+"}' mm9.fa.size.win25k.bed > mm9.fa.size.win25k.bed.SAF
awk -vOFS='\t' '{print $4,$1,$2,$3,"+"}' mm9.fa.size.win50k.bed > mm9.fa.size.win50k.bed.SAF

featureCounts -T 24 -B -C -F SAF -O -p -a  mm9.fa.size.win10k.step500.bed.SAF -o mm9.fa.size.win10k.step500.bed.ft.txt  B1-100U.mm9.sort.bam B2-150U.mm9.sort.bam B3-180U.mm9.sort.bam B4-250U.mm9.sort.bam Sr_L3_Q801601.R1_val.filter.mm9.sort.bam H3K9me3_merge.bam H3K27me3_merge.bam largeMerge.mm9.bam smallMerge.mm9.bam E14input_liuY_WT.mm9.sort.bam InputMerged.bam > log.count 2>&1 &

featureCounts -T 24 -B -C -F SAF -O -p -a  mm9.fa.size.win10k.step500.bed.SAF -o mm9.fa.size.win10k.step500.bed.ft.txt  B1-100U.mm9.sort.bam B2-150U.mm9.sort.bam B3-180U.mm9.sort.bam B4_200U.mm9.sort.bam Sr_L3_Q801601.R1_val.filter.mm9.sort.bam largeMerge.mm9.bam smallMerge.mm9.bam E14input_liuY_WT.mm9.sort.bam InputMerged.bam > log.count 2>&1 &

featureCounts -T 24 -B -C -F SAF -O -p -a  mm9.fa.size.win10k.bed.SAF -o mm9.fa.size.win10k.bed.ft.B1234DNaseMNase.txt E14input_liuY_WT.mm9.sort.bam  DNase_merge.bam B1-100U.mm9.sort.bam  B2-150U.mm9.sort.bam  B3-180U.mm9.sort.bam  B4-250U.mm9.sort.bam > log.B1234DNaseMNase.10k.count 2>&1 &


featureCounts -T 24 -B -C -F SAF -O -p -a  mm9.fa.size.win25k.bed.SAF -o mm9.fa.size.win25k.bed.ft.B1234DNaseMNase.txt E14input_liuY_WT.mm9.sort.bam  DNase_merge.bam B1-100U.mm9.sort.bam  B2-150U.mm9.sort.bam  B3-180U.mm9.sort.bam  B4-250U.mm9.sort.bam > log.B1234DNaseMNase.25k.count 2>&1 &

featureCounts -T 24 -B -C -F SAF -O -p -a  mm9.fa.size.win50k.bed.SAF -o mm9.fa.size.win50k.bed.ft.B1234DNaseMNase.txt E14input_liuY_WT.mm9.sort.bam  DNase_merge.bam B1-100U.mm9.sort.bam  B2-150U.mm9.sort.bam  B3-180U.mm9.sort.bam  B4-250U.mm9.sort.bam > log.B1234DNaseMNase.50k.count 2>&1 &



##http://bioinf.wehi.edu.au/featureCounts/
featureCounts -T 6 -B -C -F SAF -O -p -a  mm9.fa.size.win10k.step500.bed.SAF -o B4_counts.ft.txt  B4_200U.mm9.sort.bam

#-p: treat as fragment (will automaticly detect SE or PE)
#can count on all bam files SE PE mixed
#-F simple annotation file
#-O eq to Rsubread allowMultiOverlap, useful when counting overlap sliding windows

#cmd featureCount result eq Rsubread featureCount



