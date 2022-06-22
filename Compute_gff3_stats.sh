gffcompare="./../../../../../Descargas/gffcompare/gffcompare"
reference_cds="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/protein_matches/Drosophila_melanogaster.BDGP6.32.105.gff3.for_eval.over200"
sequence="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.curated"
out_folder="/home/jaume/Escritorio/Internship_CRG/training/invertebrates/prot_dmnd/output/"
rm ${out_folder}*.loci; rm ${out_folder}*.tracking; rm ${out_folder}*.annotated.gtf;
# rm output/*.stats
# ALL MATCHING PROTEINS WITHOUT LOOKING AT ORFs (full CDS sequences)
echo "Computing metrics on ${sequence}";
echo "*********************************";
echo "Prot_matching_thresh: ALL";
${gffcompare} -T -r ${reference_cds} ${sequence}.gff3 -o ${sequence}.stats
echo "*********************************";
rm ${out_folder}*.loci; rm ${out_folder}*.tracking; rm ${out_folder}*.annotated.gtf;
# MATCHING PROTEINS OVER A THRESHOLD VALUE WITHOUT LOOKING AT ORFs (full CDS sequences)
for threshold in 100 200 400 600 800; do
#     #awk '$6>800' output/Homo_sapiens.GRCh38.dna.primary_assembly.Vertebrates.5references.curated.gff3 | cut -f-9 > output/Homo_sapiens.GRCh38.dna.primary_assembly.Vertebrates.5references.curated.over800.gff3
     echo "*********************************";
     echo "Prot_matching_thresh: ALL Min_orf_len: ${min_orf_len}";
     ${gffcompare} -T -r ${reference_cds} ${sequence}.over${threshold}.gff3 -o ${sequence}.over${threshold}.stats
     echo "*********************************";
done;
rm ${out_folder}*.loci; rm ${out_folder}*.tracking; rm ${out_folder}*.annotated.gtf;
# ALL MATCHING PROTEINS LOOKING AT DIFFERENT ORF LENGTHS
for min_orf_len in 40 60 90 100 120; do
    echo "*********************************";
    echo "Prot_matching_thresh: ALL Min_orf_len: ${min_orf_len}";
    ${gffcompare} -T -r ${reference_cds} ${sequence}.ORF.${min_orf_len}_longest.gff3 -o ${sequence}.ORF.${min_orf_len}_longest.stats
    echo "*********************************";
done;
rm ${out_folder}*.loci; rm ${out_folder}*.tracking; rm ${out_folder}*.annotated.gtf;
# MATCHING PROTEINS OVER A THRESHOLD VALUE AT DIFFERENT ORF LENGTHS
for threshold in 100 200 400 600 800; do
    for min_orf_len in 40 60 90 100 120; do
        echo "*********************************";
        echo "Prot_matching_thresh: ${threshold} Min_orf_len: ${min_orf_len}";
        ${gffcompare} -T -r ${reference_cds} ${sequence}.over${threshold}.ORF.${min_orf_len}_longest.gff3 -o ${sequence}.over${threshold}.ORF.${min_orf_len}_longest.stats
        echo "*********************************";
    done;
done;
rm ${out_folder}*.loci; rm ${out_folder}*.tracking; rm ${out_folder}*.annotated.gtf;
