## Analyze tcs per genome

for efile in ./results/obtain_and_process_tcs_output/*all.tsv; do
	efile=${efile##*/}
	./pipeline/analyze_tcs_per_genome.py -i ${efile} -s ./input/MiST_domains_18.tsv -f \
	./results/analyze_tcs_per_genome_output/${efile%.*}_domains.tsv -g \
	./results/analyze_tcs_per_genome_output/${efile%.*}_domain_comb.tsv -k ./results/analyze_tcs_per_genome_output/${efile%.*}_superfamily.tsv -l ./results/analyze_tcs_per_genome_output/${efile%.*}_superfamily_comb.tsv
done


##

