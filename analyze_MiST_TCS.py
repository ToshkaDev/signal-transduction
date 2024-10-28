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
	-i || --ifile              - input file
'''

#Variables controlled by the script parameters
INPUT_FILE = None
OUTPUT_FILE = "output_HK.tsv"

#Variables set within the script
LOGGER = logging.getLogger(__name__)
logging.basicConfig(filename=sys.argv[0].replace(".py", "") + "_log.txt", level=logging.INFO)

def initialize(argv):
	global INPUT_FILE, OUTPUT_FILE
	try:
		opts, args = getopt.getopt(argv[1:],"hi:o:",["help", "ifile=", "ofile="])
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
				INPUT_FILE = str(arg).strip()
			elif opt in ("-o", "--ofile"):
				OUTPUT_FILE = str(arg).strip()
	except Exception as e:
		print("===========ERROR==========\n " + str(e) + USAGE)
		sys.exit(2)

HIS_KINASE_DIM_DOMAINS = ["HisKA", "HisKA_2", "HisKA_3", "H-kinase_dim", "His_kinase"]
HIS_KINASE_CATAL_DOMAINS = ["HATPase_c", "HATPase_c_2", "HATPase_c_5", "HWE_HK"]
RESPONSE_REG_DOMAINS = ["Response_reg", "FleQ"]


#{"GCF_xxx": {}}



#Genome_id,NCBI_id,MiST_id,protein_length,domain_architecture,sensors_or_regulators,domain_counts,domain_combinations = [f for f in protein.strip().split("\t")]


GENOME_TO_DOMAIN_DATA = {}
GENOME_TO_DOMAIN_COMB_DATA = {}


# domain_to_protein_count = {"dCache_1":5, "GAF": 2, ...}. This means that in a given genome dCache was found in 5 proteins and GAF in 2 of a considered system.
# This does not count how many times a domain is present in a single protein.
# domain_comb_to_protein_count = {"MEDS,PAS,PAS_3,PAS_4,PAS_9":2}. This means that a particualr domain combiation was found in 2 protins in the given genome
def findDomainAndCombPrevalenceInProteins():
	domain_to_protein_count = collections.defaultdict(int)
	domain_comb_to_protein_count = collections.defaultdict(int) 
	sensor_domain_comb = ""

	Genome_id_prev = None
	Genome_id = "Holder"
	with open(INPUT_FILE, "r") as iFile:
		for protein in iFile:
			#filed 6 has only uniqe domain names with indicated counts showing how many times a given domain is present in a given protein.
			# Ex., HATPase_c:1,HisKA_2:1,MEDS:1,PAS:1,PAS_3:2,PAS_4:1,PAS_9:1 
			domain_counts = protein.replace("<", "").replace(">", "").strip().split("\t")[6]
			#if records of a new genome began, save the previous genome data and update variables
			if Genome_id != Genome_id_prev:
				if Genome_id_prev:
					GENOME_TO_DOMAIN_DATA[Genome_id_prev] = domain_to_protein_count
					GENOME_TO_DOMAIN_COMB_DATA[Genome_id_prev] = domain_comb_to_protein_count
				Genome_id_prev = Genome_id
				domain_to_protein_count = collections.defaultdict(int)
				domain_comb_to_protein_count = collections.defaultdict(int)

			domains_and_counts = domain_counts.split(",")
			for domain in domains_and_counts:
				domain = domain.split(":")[0]
				#if entry is a sensor domain
				if domain not in HIS_KINASE_DIM_DOMAINS and domain not in HIS_KINASE_CATAL_DOMAINS:
					domain_to_protein_count[domain]+=1
					sensor_domain_comb = sensor_domain_comb + "," + domain
			#count this domain combination
			domain_comb_to_protein_count[sensor_domain_comb]+=1 
			

def main(argv):
	initialize(argv)
	findDomainAndCombPrevalenceInProteins()

main(sys.argv)
