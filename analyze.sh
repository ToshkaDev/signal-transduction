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
		grep "$db" ${BACT_FILE} > ${BACT_FILE%.*}_"${DB[$db]}".tsv
		grep "$db" ${ARCH_FILE} > ${ARCH_FILE%.*}_"${DB[$db]}".tsv
	done
	
	obtain_and_prepare_gtdb_files

	if [ -f "./results/obtain_and_process_tcs/output.tar.gz" ]; then
		echo "Unpack and decompress the results of querying mistdb.com to obtain domain information for analyzed genomes" 
		tar xvf ./results/obtain_and_process_tcs/output.tar.gz -C ./results/obtain_and_process_tcs/
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

		# Extract from the GTDB metadat file only those records that comrreposnd to genomes in repr_set_v214_Oct2024_MiST_MetaMiST_s.tsv file:
		# join based on the first field (genome version field) of the epr_set_v214_Oct2024_MiST_MetaMiST_s.tsv
		join -1 1 -t $'\t' ${REPR_FILE%.*}_s.tsv ${COMBINED_FILE}_${VERSION}_s.tsv > ${COMBINED_FILE}_${VERSION}.tsv

		# Extract needed fields, remove GB_ and RS_ parts from the GTDB genome ID field ($1), and add the genome_accession field (gacc[1]).
		# Fields: $14 - genome_size, #89 - protein count, #17 - GDTB taxonomy; $79 - NCBI taxonomy.
		# The resulting file will be used in the pipeline.
		awk 'BEGIN{FS="\t"; OFS="\t"} {split($1, gacc2, "_"); version=gacc2[2]"_"gacc2[3]; split(version, gacc, "."); print version,gacc[1],$14,$89,$17,$79}' ${COMBINED_FILE}_${VERSION}.tsv > ${COMBINED_FILE}_${VERSION}_p.tsv

		# Prepare a file for the database: replace ';' by tabs, remove from GTDB taxonomy 'd__', 'p__', etc. paret (see this d__Archaea;p__Nanoarchaeota;...).
		# The resuling file will be stored in the database.
		sed -e 's/;/\t/g' -e 's/[[:alpha:]]__//g' ${COMBINED_FILE}_${VERSION}_p.tsv > ${COMBINED_FILE}_${VERSION}_db.tsv

		# remove unused files 
		rm ${METADA_AR_FILE}$EXT  ${METADATA_BAC_FILE}$EXT ${METADA_AR_FILE}.tsv ${METADATA_BAC_FILE}.tsv ${COMBINED_FILE}_${VERSION}_s.tsv ${COMBINED_FILE}_${VERSION}.tsv
	fi
}

initialize_scripts_and_folders() {
	OBTAIN="obtain_and_process_tcs.py"
	OFOLDER="./results/${OBTAIN%.*}"
	check_create "${OFOLDER}"

	ANALYZEG="analyze_tcs_per_genome.py"
	AGFOLDER="./results/${ANALYZEG%.*}"
	check_create "${AGFOLDER}"

	ANALYZET="analyze_tcs_per_taxon.py"
	ATFOLDER="./results/${ANALYZET%.*}"
	check_create "${ATFOLDER}"
}

# Obtain and perform first step analysis of two-component systems (hk - histidine kinase, rr - response regulator)
obtain() {
	echo "1. Fetching two-component systems (TCS) from MiST ..."
	echo "1.1. Fetching archaeal two-component systems ..."
	for db in ${DB[@]}; do
		./pipeline/${OBTAIN} \
			-d mistdb \
			-i ${ARCH_FILE%.*}_$db.tsv \
			-f ${OFOLDER}/hk_archaea_$db.tsv \
			-s ${OFOLDER}/rr_archaea_$db.tsv \
			-b $db
		# Put results from $DB into one file for each his kinase (hk) and resp regulators (rr)
		sed '1d' ${OFOLDER}/hk_archaea_$db.tsv >> ${OFOLDER}/hk_archaea_all.tsv
		sed '1d' ${OFOLDER}/rr_archaea_$db.tsv >> ${OFOLDER}/rr_archaea_all.tsv
	done

	echo "1.2. Fetching bacterial two-component systems ..."
	for db in ${DB[@]}; do
		./pipeline/${OBTAIN} \
			-d mistdb \
			-i ${BACT_FILE%.*}_$db.tsv \
			-f ${OFOLDER}/hk_bacteria_$db.tsv \
			-s ${OFOLDER}/rr_bacteria_$db.tsv \
			-b $db
		# Put results from $DB into one file for each his kinase (hk) and resp regulators (rr)
		sed '1d' ${OFOLDER}/hk_bacteria_$db.tsv >> ${OFOLDER}/hk_bacteria_all.tsv
		sed '1d' ${OFOLDER}/rr_bacteria_$db.tsv >> ${OFOLDER}/rr_bacteria_all.tsv
	done
	cat ${OFOLDER}/hk_archaea_all.tsv ${OFOLDER}/rr_archaea_all.tsv ${OFOLDER}/hk_bacteria_all.tsv ${OFOLDER}/rr_bacteria_all.tsv > ${OFOLDER}/per_protein_combined_db.tsv
}

analyze() {
	if [ -z "$(ls -A $AGFOLDER)" ]; then
		analyze_systems_by_genome
	else
		echo "Folder '$AGFOLDER' is not empty. Analysis of two-component systems by genome will not run. Empty the folder before running the pipeline."
	fi

	if [ -z "$(ls -A $ATFOLDER)" ]; then
		analyze_systems_by_taxon
	else
		echo "Folder '$ATFOLDER' is not empty.  Analysis of two-component systems by taxon will not run. Empty the folder before running the pipeline."
	fi
}

# Analyze two-component systems by genome
analyze_systems_by_genome() {
	echo "2. Analyzing two-component systems by genome ..."
	for efile in ${OFOLDER}/*all.tsv; do
		edfile=${efile##*/}
		# One of: 'hk', 'rr', 'ocp'
		ptype=${edfile%%_*}
		./pipeline/${ANALYZEG} \
			-i ${efile} \
			-s ./input/MiST_domains_18.tsv \
			-d mistdb \
			-p $ptype \
			-t ${COMBINED_FILE}_${VERSION}_p.tsv \
			-f ${AGFOLDER}/${edfile%.*}_domains.tsv \
			-g ${AGFOLDER}/${edfile%.*}_domain_comb.tsv \
			-k ${AGFOLDER}/${edfile%.*}_superfamily.tsv \
			-l ${AGFOLDER}/${edfile%.*}_superfamily_comb.tsv
	done
	cat ${AGFOLDER}/*.tsv > ${AGFOLDER}/per_genome_combined_db.tsv
}

# Analyze two-component systems by taxonomic level
analyze_systems_by_taxon() {
	echo "3. Analyzing two-component systems by taxa using files generated in the previous step ..."
	levels=("species" "genus" "family" "order" "class" "phylum" "kingdom")
	for level in ${levels[@]}; do
		for efile in ${AGFOLDER}/*.tsv; do
			# ./results/hk_archaea_all_superfamily_comb.tsv -> hk_archaea_all_superfamily_comb.tsv
			edfile=${efile##*/}
			# hk_archaea_all_superfamily_comb.tsv -> hk
			ptype=${edfile%%_*}
			# hk_archaea_all_superfamily_comb.tsv -> hk_archaea_all_superfamily_comb
			subname=${edfile%.*}
			# hk_archaea_all_superfamily_comb -> superfamily_comb
			ctype=${subname#*all}
			./pipeline/${ANALYZET} \
				-i ${efile} \
				-s ${COMBINED_FILE}_${VERSION}_p.tsv \
				-d mistdb \
				-p $ptype \
				-c $ctype \
				-f ${ATFOLDER}/${edfile%.*}_$level.tsv \
				-t $level
		done
	done
	cat ${ATFOLDER}/*.tsv > ${ATFOLDER}/per_taxon_combined_db.tsv
}

prepare_files
initialize_scripts_and_folders
# Run the process 'obtain' only if $OFOLDER is empty
[ -z "$(ls -A $OFOLDER)" ] && obtain
analyze
