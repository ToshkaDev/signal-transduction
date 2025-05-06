#!/bin/bash

declare -A DB=(["MiST"]="mist" ["MetaMiST"]="mist-mags")

check_create() {
	[ ! -d $1 ] && mkdir $1
}

prepare_files() {
	BACT_FILE="./input/repr_set_v214_Oct2024_MiST_MetaMiST_bact.tsv"
	ARCH_FILE="./input/repr_set_v214_Oct2024_MiST_MetaMiST_arch.tsv"

	grep "d__Bacteria" ${BACT_FILE%_*}.tsv > ${BACT_FILE}
	grep "d__Archaea" ${ARCH_FILE%_*}.tsv > ${ARCH_FILE}

	for db in ${!DB[@]}; do
		grep "$db" ${BACT_FILE} > ${BACT_FILE%.*}_"${DB[$db]}".tsv
		grep "$db" ${ARCH_FILE} > ${ARCH_FILE%.*}_"${DB[$db]}".tsv
	done

	echo "Decompress GTDB metadata file" 
	unzip ./input/gtdb_metadata/ar53_bac120_metadata_r214.tsv.zip -d ./input/gtdb_metadata
	
	echo "Unpack and decompress the results of querying mistdb.com to obtain domain information for analyzed genomes" 
	tar xvf ./results/obtain_and_process_tcs/output.tar.gz -C ./results/obtain_and_process_tcs/
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

# Obtain and perform first step analysis of two-component systems
obtain() {
	echo "Starting 'obtain' ..."
	echo "Fetching two-component systems (TCS) from MiST ..."
	echo "Fetching archaeal ..."
	for db in ${DB[@]}; do
		./pipeline/${OBTAIN} -i ${ARCH_FILE%.*}_$db.tsv \
			-f ${OFOLDER}/his_kinases_archaea_$db.tsv \
			-s ${OFOLDER}/resp_regulators_archaea_$db.tsv \
			-d $db
		# put results from $DB into one file for each his_kinase and resp_regulators
		cat ${OFOLDER}/his_kinases_archaea_$db.tsv >> ${OFOLDER}/his_kinases_archaea_all.tsv
		cat ${OFOLDER}/resp_regulators_archaea_$db.tsv >> ${OFOLDER}/resp_regulators_archaea_all.tsv
	done

	echo "Fetching bacterial two-component systems ..."
	for db in ${DB[@]}; do
		./pipeline/${OBTAIN} -i ${BACT_FILE%.*}_$db.tsv \
			-f ${OFOLDER}/his_kinases_bacteria_$db.tsv \
			-s ${OFOLDER}/resp_regulators_bacteria_$db.tsv \
			-d $db
		# put results from $DB into one file for each his_kinase and resp_regulators
		cat ${OFOLDER}/his_kinases_bacteria_$db.tsv >> ${OFOLDER}/his_kinases_bacteria_all.tsv
		cat ${OFOLDER}/resp_regulators_bacteria_$db.tsv >> ${OFOLDER}/resp_regulators_bacteria_all.tsv
	done
}

# Analyze two-component systems per genome and per taxonnomic level
analyze() {
	echo "Starting 'analyze' ..."
	echo "Analyzing two-component systems per genome ..."
	for efile in ${OFOLDER}/*all.tsv; do
		edfile=${efile##*/}
		./pipeline/${ANALYZEG} -i ${efile} -s ./input/MiST_domains_18.tsv \
		-t ./input/gtdb_metadata/ar53_bac120_metadata_r214_ed.tsv \
		-f ${AGFOLDER}/${edfile%.*}_domains.tsv \
		-g ${AGFOLDER}/${edfile%.*}_domain_comb.tsv \
		-k ${AGFOLDER}/${edfile%.*}_superfamily.tsv \
		-l ${AGFOLDER}/${edfile%.*}_superfamily_comb.tsv
	done

	echo "Analyzing two-component systems per taxon using the files generated at the previous step ..."
	levels=("species" "genus" "family" "order" "class" "phylum" "kingdom")
	for level in ${levels[@]}; do
		for efile in ${AGFOLDER}/*.tsv; do
			edfile=${efile##*/}
			./pipeline/${ANALYZET} -i ${efile} -s ./input/gtdb_metadata/ar53_bac120_metadata_r214_ed.tsv \
			-f ${ATFOLDER}/${edfile%.*}_$level.tsv \
			-t $level
		done
	done	
}

prepare_files
initialize_scripts_and_folders
# run the process 'obtain' only if $OFOLDER is empty
[ -z "$(ls -A $OFOLDER)" ] && obtain
analyze


