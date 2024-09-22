
cat select.up.log2TPM.ESC.mat|perl mat_add_desc.pl refSeq.ucsc.name.uniq.description > select.up.log2TPM.ESC.mat.txt
diff -s <(cut -f 1 select.up.log2TPM.ESC.mat.txt) <(cut -f 1 select.up.log2TPM.ESC.mat) |les
#identifal

fisher.pl cluster3/expected.up.list select.up.log2TPM.ESC.mat


#cut -f 1 select.up.log2TPM.ESC.mat|sed '/^$/d; s/\s+//g'|sort > select.up.log2TPM.ESC.mat.id

#fisher.pl select.up.log2TPM.ESC.mat.id refSeq.ucsc.name.uniq.description > select.up.log2TPM.ESC.mat.id.decs

#diff -s <(sort select.up.log2TPM.ESC.mat.id) <(cut -f 1 select.up.log2TPM.ESC.mat.id.decs|sort ) |les #identical



