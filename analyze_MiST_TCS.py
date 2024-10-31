#!/usr/bin/python3
import sys, getopt
import urllib.request, urllib.parse, urllib.error
import json
import collections
import os.path
import time
import logging

USAGE = "\n\nThe script . \n" + \
	". \n" + \
	"python 	" + sys.argv[0] + '''
	-h || --help               - help
	-i || --ifile              - input file 1
	-s || --sfile              - input file 2 (with MiST domain information)
	-f || --ffile              - output file 1
	-g || --gfile              - output file 2 
'''

# Variables controlled by the script parameters
INPUT_FILE1 = None
INPUT_FILE2 = None
OUTPUT_FILE1 = "genome_to_domain_data.tsv"
OUTPUT_FILE2 = "genome_to_domain_comb_data.tsv"

MIST_DOMAIN_TO_SUPERFAMILY = {}
HIS_KINASE_DIM_DOMAINS = ["HisKA", "HisKA_2", "HisKA_3", "H-kinase_dim", "His_kinase"]
HIS_KINASE_CATAL_DOMAINS = ["HATPase_c", "HATPase_c_2", "HATPase_c_5", "HWE_HK"]
RESPONSE_REG_DOMAINS = ["Response_reg", "FleQ"]

# Variables set within the script
LOGGER = logging.getLogger(__name__)
logging.basicConfig(filename=sys.argv[0].replace(".py", "") + "_log.txt", level=logging.INFO)

def initialize(argv):
	global INPUT_FILE1, INPUT_FILE2, OUTPUT_FILE1, OUTPUT_FILE2
	try:
		opts, args = getopt.getopt(argv[1:],"hi:s:f:g:",["help", "ifile=", "sfile=", "ffile=", "gfile="])
		if len(opts) == 0:
			raise getopt.GetoptError("Options are required\n")
	except getopt.GetoptError as e:
		print("===========ERROR==========\n " + str(e) + USAGE)
		sys.exit(2)
	try:
		for opt, arg in opts:
			if opt in ("-h", "--help"):
				print(USAGE)
				sys.exit()
			elif opt in ("-i", "--ifile"):
				INPUT_FILE1 = str(arg).strip()
			elif opt in ("-s", "--sfile"):
				INPUT_FILE1 = str(arg).strip()
			elif opt in ("-f", "--ffile"):
				OUTPUT_FILE1 = str(arg).strip()
			elif opt in ("-g", "--gfile"):
				OUTPUT_FILE2 = str(arg).strip()

	except Exception as e:
		print("===========ERROR==========\n " + str(e) + USAGE)
		sys.exit(2)

def processInput():
	global MIST_DOMAIN_TO_SUPERFAMILY
	with open(INPUT_FILE2, "r") as iFile2:
		for line in iFile2:
			domain_info = line.split("\t")
			MIST_DOMAIN_TO_SUPERFAMILY[domain_info[0]] = domain_info[3] 

#Genome_id,NCBI_id,MiST_id,protein_length,domain_architecture,sensors_or_regulators,domain_counts,domain_combinations = [f for f in protein.strip().split("\t")]

# domain_to_protein_count = {"dCache_1":5, "GAF": 2, ...}. This means that in a given genome dCache was found in 5 proteins and GAF in 2 of a considered system.
# This does not count how many times a domain is present in a single protein.
# domain_comb_to_protein_count = {"MEDS,PAS,PAS_3,PAS_4,PAS_9":2}. This means that a particualr domain combiation was found in 2 protins in the given genome
def findDomainAndCombPrevalenceInProteins():
	domain_to_protein_count = collections.defaultdict(int)
	superfamily_to_protein_count = collections.defaultdict(int)

	domain_comb_to_protein_count = collections.defaultdict(int)
	superfamily_comb_to_protein_count = collections.defaultdict(int)

	sensor_domain_comb = ""
	sensor_superfamily_comb = {}
	Genome_id_prev = None
	with open(INPUT_FILE1, "r") as iFile:
		for protein in iFile:
			# filed 6 has only uniqe domain names with indicated counts showing how many times a given domain is present in a given protein.
			# Domains are sorted alphabetically
			# Ex., HATPase_c:1,HisKA_2:1,MEDS:1,PAS:1,PAS_3:2,PAS_4:1,PAS_9:1 
			protein_record = protein.strip().split("\t")
			Genome_id = protein_record[0]
			domain_counts = protein_record[6].replace("<", "").replace(">", "")
			sensor_domain_comb = ""
			# if records of a new genome began, save the previous genome data and update variables
			if Genome_id != Genome_id_prev:
				if Genome_id_prev:
					with open (OUTPUT_FILE1, "a") as oFile1, open(OUTPUT_FILE2, "a") as oFile2:
						for domain,count in domain_to_protein_count.items():
							oFile1.write("\t".join([Genome_id_prev, domain, str(count)]) + "\n")
						for domain_comb,count in domain_comb_to_protein_count.items():
							oFile2.write("\t".join([Genome_id_prev, domain_comb, str(count)]) + "\n")

				Genome_id_prev = Genome_id
				domain_to_protein_count = collections.defaultdict(int)
				domain_comb_to_protein_count = collections.defaultdict(int)

			domains_and_counts = domain_counts.split(",")
			for domain in domains_and_counts:
				domain = domain.split(":")[0]
				# if entry is a sensor domain
				if domain not in HIS_KINASE_DIM_DOMAINS and domain not in HIS_KINASE_CATAL_DOMAINS:
					domain_to_protein_count[domain]+=1
					superfamily_to_protein_count[MIST_DOMAIN_TO_SUPERFAMILY[domain]]+=1

					sensor_domain_comb = sensor_domain_comb + "," + domain
					sensor_superfamily_comb[MIST_DOMAIN_TO_SUPERFAMILY[domain]] = None
			# count this domain combination
			domain_comb_to_protein_count[sensor_domain_comb.lstrip(",")]+=1 
			

def main(argv):
	initialize(argv)
	findDomainAndCombPrevalenceInProteins()

main(sys.argv)
