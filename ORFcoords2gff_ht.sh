sequence="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.curated"
reference_fa="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.fa"
out_folder="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/"
echo "Add the coords in the curated gff3"
for threshold in 100 200 400 600 800; do
	for min_orf_len in 40 60 90 100 120; do
		echo "Threshold: ${threshold} ORF: ${min_orf_len}"
		python3 ORFcoords2gff.py ${sequence}.over${threshold}.getorf.${min_orf_len}_longest.bed ${sequence}.over${threshold}.gff3 ${sequence}.over${threshold}.ORF.${min_orf_len}_longest.gff3
	done;
done;
echo "****************************************************"
for min_orf_len in 40 60 90 100 120; do
	echo "ORF: ${min_orf_len}"
	python3 ORFcoords2gff.py ${sequence}.getorf.${min_orf_len}_longest.bed ${sequence}.gff3 ${sequence}.ORF.${min_orf_len}_longest.gff3
done;
echo "****************************************************"
