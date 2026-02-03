#!/bin/bash

declare -A DB=(["MiST"]="mist" ["MetaMiST"]="mist-mags")
REPR_FILE="./input/repr_set_v214_Oct2024_MiST_MetaMiST.tsv"

check_create() {
	[ ! -d $1 ] && mkdir $1
}

prepare_files() {
	echo "Started downloading and preparing files ..."
	BACT_FILE=${REPR_FILE%.*}_bact.tsv
	ARCH_FILE=${REPR_FILE%.*}_arch.tsv

	grep "d__Bacteria" ${REPR_FILE} > ${BACT_FILE}
	grep "d__Archaea" ${REPR_FILE} > ${ARCH_FILE}

	for db in ${!DB[@]}; do
		grep "\s$db" ${BACT_FILE} > ${BACT_FILE%.*}_"${DB[$db]}".tsv
		grep "\s$db" ${ARCH_FILE} > ${ARCH_FILE%.*}_"${DB[$db]}".tsv
	done
	
	obtain_and_prepare_gtdb_files

	OFOLDER="./results/obtain_and_process_st"
	if [ -f "${OFOLDER}/tcs.tar.gz" ] || [ -f "${OFOLDER}/ocp_bacteria_all1.tsv.tar.gz" ]; then
		echo "Unpack and decompress the results of querying mistdb.com to obtain domain information for analyzed genomes" 
		tar xvf ${OFOLDER}/tcs.tar.gz -C ${OFOLDER}/
		tar xvf ${OFOLDER}/ocp_archaea_all.tsv.tar.gz -C ${OFOLDER}/
		tar xvf ${OFOLDER}/ocp_bacteria_all1.tsv.tar.gz -C ${OFOLDER}/
		tar xvf ${OFOLDER}/ocp_bacteria_all2.tsv.tar.gz -C ${OFOLDER}/
		tar xvf ${OFOLDER}/ocp_bacteria_all3.tsv.tar.gz -C ${OFOLDER}/
        tar xvf ${OFOLDER}/mcp.tar.gz -C ${OFOLDER}/
		# Concatenate ocp_bacteria_all* files
		cat ${OFOLDER}/ocp_bacteria_all1.tsv ${OFOLDER}/ocp_bacteria_all2.tsv ${OFOLDER}/ocp_bacteria_all3.tsv > ${OFOLDER}/ocp_bacteria_all.tsv

		# Prepare the database file
		echo $'genome\tgenome_accession\tncbi_protein_accession\tmist_protein_accession\tprotein_type\tsource\tprotein_length\t'\
			$'domain_architecture\tsensors_or_regulators\tdomain_counts\tdomains' | sed 's/ //g' > ${OFOLDER}/per_protein_combined_db.tsv
		cat ${OFOLDER}/hk_archaea_all.tsv ${OFOLDER}/rr_archaea_all.tsv ${OFOLDER}/hk_bacteria_all.tsv ${OFOLDER}/rr_bacteria_all.tsv \
			${OFOLDER}/ocp_archaea_all.tsv ${OFOLDER}/ocp_bacteria_all.tsv \
            ${OFOLDER}/mcp_archaea_all.tsv ${OFOLDER}/mcp_bacteria_all.tsv >> ${OFOLDER}/per_protein_combined_db.tsv
	fi
}

obtain_and_prepare_gtdb_files() {
	VERSION=r214
	COMBINED_FILE=./input/gtdb_metadata/ar_bac_metadata

	if [ ! -f "${COMBINED_FILE}_${VERSION}_p.tsv" ] || [ ! -f "${COMBINED_FILE}_${VERSION}_db.tsv" ]; then
		METADA_AR_LINK=https://data.gtdb.ecogenomic.org/releases/release214/214.0/ar53_metadata_r214.tar.gz
		METADATA_BAC_LINK=https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz
		METADA_AR_FILE=./input/gtdb_metadata/ar_metadata
		METADATA_BAC_FILE=./input/gtdb_metadata/bac_metadata
		EXT=.tar.gz

		echo "Downloading the Genome Taxonomy Database (version ${VERSION}) metadata files ..."
		wget -O ${METADA_AR_FILE}$EXT ${METADA_AR_LINK}
		wget -O ${METADATA_BAC_FILE}$EXT ${METADATA_BAC_LINK}

		echo "Decompressing metadata files ..."
		tar xOvf ${METADA_AR_FILE}$EXT | sed '1d' > ${METADA_AR_FILE}.tsv
		tar xOvf ${METADATA_BAC_FILE}$EXT | sed '1d' > ${METADATA_BAC_FILE}.tsv

		echo "Preparing metadata files ..."
		# Combine both files
		cat ${METADA_AR_FILE}.tsv ${METADATA_BAC_FILE}.tsv > ${COMBINED_FILE}_${VERSION}.tsv

		# Sort genomic files to prepare for use with the 'join' operator below
		cut -f 1 ${REPR_FILE} | sort  > ${REPR_FILE%.*}_s.tsv
		sort -k 1 ${COMBINED_FILE}_${VERSION}.tsv > ${COMBINED_FILE}_${VERSION}_s.tsv

		# Extract from the GTDB metadat file only those records that correspond to genomes in repr_set_v214_Oct2024_MiST_MetaMiST_s.tsv file:
		# join based on the first field (genome version field) of the epr_set_v214_Oct2024_MiST_MetaMiST_s.tsv
		join -1 1 -t $'\t' ${REPR_FILE%.*}_s.tsv ${COMBINED_FILE}_${VERSION}_s.tsv > ${COMBINED_FILE}_${VERSION}.tsv

		# Extract needed fields, remove GB_ and RS_ parts from the GTDB genome ID field ($1), and add the genome_accession field (gacc[1]).
		# Fields: $14 - genome_size, #89 - protein count, #17 - GDTB taxonomy; $79 - NCBI taxonomy.
		# The resulting file will be used in the pipeline.
		awk 'BEGIN{FS="\t"; OFS="\t"} {split($1, gacc2, "_"); version=gacc2[2]"_"gacc2[3]; split(version, gacc, "."); print version,gacc[1],$14,$89,$17,$79}' ${COMBINED_FILE}_${VERSION}.tsv > ${COMBINED_FILE}_${VERSION}_p.tsv

		# Prepare a file for the database: add headers, replace ';' by tabs, remove from GTDB taxonomy 'd__', 'p__', etc. paret (see this d__Archaea;p__Nanoarchaeota;...).
		# sed 's/ //g' is needed because echo with the multiline string introduces a space with the string continuation at the new line (as a result of string formating by VS)
		# The resuling file will be stored in the database.
		echo $'genome_version\tgenome_accession\tgenome_size\tprotein_count\tgtdb_kingdom\tgtdb_phylum\tgtdb_class\tgtdb_order\tgtdb_family\t'\
			$'gtdb_genus\tgtdb_species\tncbi_kingdom\tncbi_phylum\tncbi_class\tncbi_order\tncbi_family\tncbi_genus\tncbi_species' | sed 's/ //g' > ${COMBINED_FILE}_${VERSION}_db.tsv
		sed -e 's/;/\t/g' -e 's/[[:alpha:]]__//g' ${COMBINED_FILE}_${VERSION}_p.tsv >> ${COMBINED_FILE}_${VERSION}_db.tsv
		
		# remove unused files 
		rm ${METADA_AR_FILE}$EXT  ${METADATA_BAC_FILE}$EXT ${METADA_AR_FILE}.tsv ${METADATA_BAC_FILE}.tsv ${COMBINED_FILE}_${VERSION}_s.tsv ${COMBINED_FILE}_${VERSION}.tsv
	fi
}

initialize_scripts_and_folders() {
	OBTAIN="obtain_and_process_st.py"
	OFOLDER="./results/${OBTAIN%.*}"
	check_create "${OFOLDER}"

	ANALYZEG="analyze_st_per_genome.py"
	AGFOLDER="./results/${ANALYZEG%.*}"
	check_create "${AGFOLDER}"

	ANALYZET="analyze_st_per_taxon.py"
	ATFOLDER="./results/${ANALYZET%.*}"
	check_create "${ATFOLDER}"
}

# Obtain and perform the first step analysis of:
# two-component (hk - histidine kinase, rr - response regulator), one-component, and chemotaxis systems (mcp - methyl-accepting chemotaxis proteins)
obtain() {
	echo "1. Fetching signal transduction systems (ST) from MiST ..."
	echo "1.1. Fetching archaeal signal transduction systems ..."
	for db in ${DB[@]}; do
		./pipeline/${OBTAIN} \
			-d mistdb \
			-i ${ARCH_FILE%.*}_$db.tsv \
			-f ${OFOLDER}/hk_archaea_$db.tsv \
			-s ${OFOLDER}/rr_archaea_$db.tsv \
			-t ${OFOLDER}/ocp_archaea_$db.tsv \
            -m ${OFOLDER}/mcp_archaea_$db.tsv \
			-b $db
		# Put results from $DB into one file for each his kinase (hk) and resp regulator (rr)
		sed '1d' ${OFOLDER}/hk_archaea_$db.tsv >> ${OFOLDER}/hk_archaea_all.tsv
		sed '1d' ${OFOLDER}/rr_archaea_$db.tsv >> ${OFOLDER}/rr_archaea_all.tsv
		sed '1d' ${OFOLDER}/ocp_archaea_$db.tsv >> ${OFOLDER}/ocp_archaea_all.tsv
        sed '1d' ${OFOLDER}/mcp_archaea_$db.tsv >> ${OFOLDER}/mcp_archaea_all.tsv
	done

	echo "1.2. Fetching bacterial signal transduction systems ..."
	for db in ${DB[@]}; do
		./pipeline/${OBTAIN} \
			-d mistdb \
			-i ${BACT_FILE%.*}_$db.tsv \
			-f ${OFOLDER}/hk_bacteria_$db.tsv \
			-s ${OFOLDER}/rr_bacteria_$db.tsv \
			-t ${OFOLDER}/ocp_bacteria_$db.tsv \
            -m ${OFOLDER}/mcp_bacteria_$db.tsv \
			-b $db
		# Put results from $DB into one file for each his kinase (hk) and resp regulator (rr)
		sed '1d' ${OFOLDER}/hk_bacteria_$db.tsv >> ${OFOLDER}/hk_bacteria_all.tsv
		sed '1d' ${OFOLDER}/rr_bacteria_$db.tsv >> ${OFOLDER}/rr_bacteria_all.tsv
		sed '1d' ${OFOLDER}/ocp_bacteria_$db.tsv >> ${OFOLDER}/ocp_bacteria_all.tsv
        sed '1d' ${OFOLDER}/mcp_bacteria_$db.tsv >> ${OFOLDER}/mcp_bacteria_all.tsv
	done

	# Prepare the database file
	echo $'genome\tgenome_accession\tncbi_protein_accession\tmist_protein_accession\tprotein_type\tsource\tprotein_length\t'\
		$'domain_architecture\tsensors_or_regulators\tdomain_counts\tdomains' | sed 's/ //g' > ${OFOLDER}/per_protein_combined_db.tsv
	cat ${OFOLDER}/hk_archaea_all.tsv ${OFOLDER}/rr_archaea_all.tsv ${OFOLDER}/hk_bacteria_all.tsv ${OFOLDER}/rr_bacteria_all.tsv \
		${OFOLDER}/ocp_archaea_all.tsv ${OFOLDER}/ocp_bacteria_all.tsv ${OFOLDER}/mcp_archaea_all.tsv ${OFOLDER}/mcp_bacteria_all.tsv >> ${OFOLDER}/per_protein_combined_db.tsv

}

analyze() {
	if [ -z "$(ls -A $AGFOLDER)" ]; then
		analyze_systems_by_genome
	else
		echo "Folder '$AGFOLDER' is not empty. Analysis of signal transduction systems by genome will not run. Empty the folder before running the pipeline."
	fi

	if [ -z "$(ls -A $ATFOLDER)" ]; then
		analyze_systems_by_taxon
	else
		echo "Folder '$ATFOLDER' is not empty.  Analysis of signal transduction systems by taxon will not run. Empty the folder before running the pipeline."
	fi
}

# Analyze signal transduction systems by genome
analyze_systems_by_genome() {
	echo "2. Analyzing signal transduction systems by genome ..."
	for efile in ${OFOLDER}/*all.tsv; do
		edfile=${efile##*/}
		# One of: 'hk', 'rr', 'ocp', 'mcp'
		ptype=${edfile%%_*}
		./pipeline/${ANALYZEG} \
			-i ${efile} \
			-s ./input/MiST_domains_18.tsv \
			-d mistdb \
			-p $ptype \
			-t ${COMBINED_FILE}_${VERSION}_p.tsv \
			-f ${AGFOLDER}/${edfile%.*}_domain_p.tsv \
			-g ${AGFOLDER}/${edfile%.*}_domain_comb_p.tsv \
			-k ${AGFOLDER}/${edfile%.*}_superfamily_p.tsv \
			-l ${AGFOLDER}/${edfile%.*}_superfamily_comb_p.tsv
	done

	# Prepare the database file
	echo $'genome\tgenome_accession\tsource\tprotein_type\tdomains\tdomain_combination_type\tcount_raw\tcount_normalized_by_genome_size\t'\
		$'count_normalized_by_total_proteins' | sed 's/ //g' > ${AGFOLDER}/per_genome_combined_db.tsv
	cat ${AGFOLDER}/*p.tsv >> ${AGFOLDER}/per_genome_combined_db.tsv
}

# Analyze signal transduction systems by taxonomic level
analyze_systems_by_taxon() {
	echo "3. Analyzing signal transduction systems by taxa using files generated in the previous step ..."
	levels=("species" "genus" "family" "order" "class" "phylum" "kingdom")
	for level in ${levels[@]}; do
		for efile in ${AGFOLDER}/*p.tsv; do
			# Ex., ./results/hk_archaea_all_superfamily_comb_p.tsv -> hk_archaea_all_superfamily_comb_p.tsv
			edfile=${efile##*/}
			# hk_archaea_all_superfamily_comb_p.tsv -> hk_archaea_all_superfamily_comb.tsv
			edfile=${edfile%_*}.tsv
			# hk_archaea_all_superfamily_comb.tsv -> hk
			ptype=${edfile%%_*}
			# hk_archaea_all_superfamily_comb.tsv -> hk_archaea_all_superfamily_comb
			subname=${edfile%.*}
			# hk_archaea_all_superfamily_comb -> superfamily_comb
			ctype=${subname#*all_}
			./pipeline/${ANALYZET} \
				-i ${efile} \
				-s ${COMBINED_FILE}_${VERSION}_p.tsv \
				-d mistdb \
				-p $ptype \
				-c $ctype \
				-f ${ATFOLDER}/${edfile%.*}_${level}_p.tsv \
				-t $level
		done
	done

	# Prepare the database file
	echo $'gtdb_taxonomy_string\tgtdb_taxonomy_last\tgtdb_taxonomy_rank\tsource\tprotein_type\tdomains\tdomain_combination_type\tcount_raw\t'\
		$'count_normalized_by_total_genomes\tcount_normalized_by_genome_size_by_total_genomes\t'\
		$'count_normalized_by_total_proteins_by_total_genomes' | sed 's/ //g' > ${ATFOLDER}/per_taxon_combined_db.tsv
	cat ${ATFOLDER}/*p.tsv >> ${ATFOLDER}/per_taxon_combined_db.tsv
}

prepare_files
initialize_scripts_and_folders
# Run the process 'obtain' only if $OFOLDER is empty
[ -z "$(ls -A $OFOLDER)" ] && obtain
analyze
