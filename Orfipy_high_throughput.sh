sequence="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.curated"
reference_fa="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.fa"
out_folder="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/"

gffread -w ${sequence}.fa -g ${reference_fa} ${sequence}.gff3

for threshold in 100 200 400 600 800; do
	gffread -w ${sequence}.over${threshold}.fa -g ${reference_fa} ${sequence}.over${threshold}.gff3
done;

for min_orf_len in 40 60 90 100 120; do
    echo "*********************************";
    echo "Prot_matching_thresh: ALL Min_orf_len: ${min_orf_len}";
    orfipy ${sequence}.fa --strand f --bed ${sequence}.getorf.${min_orf_len}.bed --longest --min ${min_orf_len} --between-stops --outdir ${out_folder}
    echo "*********************************";
done;
rm ${out_folder}*.log; rm ${out_folder}*0.bed

for threshold in 100 200 400 600 800; do
    for min_orf_len in 40 60 90 100 120; do
            echo "*********************************";
            echo "Prot_matching_thresh: ${threshold} Min_orf_len: ${min_orf_len}";
            orfipy ${sequence}.over${threshold}.fa --strand f --bed ${sequence}.over${threshold}.getorf.${min_orf_len}.bed --longest --min ${min_orf_len} --between-stops --outdir ${out_folder}
            echo "*********************************";
    done;
done;

rm ${out_folder}*.log; rm ${out_folder}*0.bed
