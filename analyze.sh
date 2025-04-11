#!/bin/bash

declare -A DB=(["MiST"]="mist" ["MetaMiST"]="mist-mags")
BACT_FILE="./input/repr_set_v214_Oct2024_MiST_MetaMiST_bact.tsv"
ARCH_FILE="./input/repr_set_v214_Oct2024_MiST_MetaMiST_arch.tsv"

grep "d__Bacteria" ${BACT_FILE%_*} > ${BACT_FILE}
grep "d__Archaea" ${ARCH_FILE%_*} > ${ARCH_FILE}

for db in ${!DB[@]}; do
	grep "$db" ${BACT_FILE} > ${BACT_FILE%.*}_"$db".tsv
	grep "$db" ${ARCH_FILE} > ${ARCH_FILE%.*}_"$db".tsv
done

# Fetch TCS for Archaea:
for db in ${!DB[@]}; do
	./signal-transduction/obtain_and_process_tcs.py -i ${ARCH_FILE%.*}_"$db".tsv \
		-f ./results/obtain_and_process_tcs_output/his_kinases_archaea_"$db".tsv \
		-s ./results/obtain_and_process_tcs_output/resp_regulators_archaea_"$db".tsv \
		-d ${DB[$db]}
done

# Fetch TCS for Bacteria:
for db in ${!DB[@]}; do
	./signal-transduction/obtain_and_process_tcs.py -i ${BACT_FILE%.*}_"$db".tsv \
		-f ./results/obtain_and_process_tcs_output/his_kinases_bacteria_"$db".tsv \
		-s ./results/obtain_and_process_tcs_output/resp_regulators_bacteria_"$db".tsv \
		-d ${DB[$db]}
done


check_create() {
	[ ! -d $1 ] && mkdir $1
}

## Analyze tcs per genome
check_create "./results/analyze_tcs_per_genome_output"
for efile in ./results/obtain_and_process_tcs_output/*all.tsv; do
	edfile=${efile##*/}
	./pipeline/analyze_tcs_per_genome.py -i ${efile} -s ./input/MiST_domains_18.tsv \
	-f ./results/analyze_tcs_per_genome_output/${edfile%.*}_domains.tsv \
	-g ./results/analyze_tcs_per_genome_output/${edfile%.*}_domain_comb.tsv \
	-k ./results/analyze_tcs_per_genome_output/${edfile%.*}_superfamily.tsv \
	-l ./results/analyze_tcs_per_genome_output/${edfile%.*}_superfamily_comb.tsv
done

## Analyze tcs per taxon using the files genrated at the previous step
check_create "./results/analyze_tcs_per_taxon_output"
levels=("species" "genus" "family" "order" "class" "phylum" "kingdom")
for level in ${levels[@]}; do
	for efile in ./results/analyze_tcs_per_genome_output/*.tsv; do
		edfile=${efile##*/}
		./pipeline/analyze_tcs_per_taxon.py -t $level -i ${efile} -s ./input/gtdb_taxonomy/*.tsv \
		-f ./results/analyze_tcs_per_taxon_output/${edfile%.*}_$level.tsv
	done
done



