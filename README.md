# Geneidx-optimize-params
Second part of master's thesis focused on updating the geneid program by both introducing protein matching and optimizing the parameters file

#										**************************
#										*GENE PREDICTION PIPELINE*
#										**************************
# Preparation
	# 1. To obtain the set of proteins (insecta example):
		#/nfs/users/rg/jreig/coding_matches
		#wget "https://www.uniprot.org/uniref/?query=taxonomy:50557+AND+count:[15 TO *]+AND+identity:0.9&format=fasta&compress=yes" -O UniRef90.50557.15+.fa.gz
		#gunzip -c UniRef90.50557.15+.fa.gz > UniRef90.50557.15+.fa;
		#diamond makedb --in UniRef90.50557.15+.fa -d UniRef90.50557.15+;
	# 2. To obtain species' genome fasta (Drosophila melanogaster example):
		#/nfs/users/rg/jreig/coding_matches/data
		# website: http://ftp.ensembl.org/pub/release-105/fasta/drosophila_melanogaster/dna/
		# wget http://ftp.ensembl.org/pub/release-105/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.fa.gz
# 1) run the genome of species against a set of curated proteins

# 2) Find ORF and compute stats

# 3) Creating parameter file (auto-training)

# 4) Run nextflow geneidx pipeline
