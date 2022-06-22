#
# Author: Ferriol Calvet (ferriol.calvet@crg.eu)
#

file='Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.UniRef90.50557.15+.hsp'
path='/home/jaume/Escritorio/Internship_CRG/training/invertebrates/build_param_droso_600_90/'

sort -k1,1 -k4,5n -k9,9 ${path}${file} | awk '!found[$1"\t"$2"\t"$3"\t"$4]++' | awk '!found[$1"\t"$2"\t"$3"\t"$5]++' > ${path}${file}.summarized_matches.gff

filee="Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.hsp.summarized_matches.gff"

awk 'OFS="\t"{print $1, $4-40, $5+40, $9$7}' $path$filee | sort -k1,1 -k4,4 -k2,2n > $path$filee.resorted

python compute_introns.py $path$filee.resorted $path$filee.introns 10000

awk '!found[$1"\t"$2"\t"$3"\t"$4"\t"$5]++' $path$filee.introns | sort -k1,1 -k4,5n > $path$filee.introns.non_redundant

bedtools intersect -a ${path}Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.hsp.summarized_matches.gff.introns.non_redundant -b ${path}Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.hsp.gff -v > Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.hsp.summarized_matches.gff.introns.non_redundant.non_overlapping_matches



