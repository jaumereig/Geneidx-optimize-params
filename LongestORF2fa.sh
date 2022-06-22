sequence="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.curated"
reference_fa="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.fa"
out_folder="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/"

echo "Extract the sequence of ORF regions in fasta format to train geneid"

for threshold in 100 200 400 600 800; do
	for min_orf_len in 40 60 90 100 120; do
		echo "Threshold: ${threshold}. ORF: ${min_orf_len}"
		gffread -w ${sequence}.over${threshold}.ORF.${min_orf_len}_longest.fa -g ${reference_fa} ${sequence}.over${threshold}.ORF.${min_orf_len}_longest.gff3
	done;
done;
echo "****************************************************"
for min_orf_len in 40 60 90 100 120; do
        echo "Threshold: All. ORF: ${min_orf_len}"
	gffread -w ${sequence}.ORF.${min_orf_len}_longest.fa -g ${reference_fa} ${sequence}.ORF.${min_orf_len}_longest.gff3
done;
