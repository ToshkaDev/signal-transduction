#!/usr/bin/python3

import sys
import getopt
import urllib.request, urllib.error
import json
from collections import OrderedDict
import time
import logging

USAGE = "\n\nThe script queries MiST database via it's API for the counts of signal transduction proteins in genomes.\n\n" + \
	"It produces ST info for each genome. \n\n" + \
	"python 	" + sys.argv[0] + '''
	-h || --help               - help
	-i || --ifile              - input file
	-s || --sfile              - input file 2 (GTDB taxonomy metadata file)
	-o || --ofile              - output file with major st systems (one-component, two-component, ect, other)
	-c || --cfile			   - output file with chemotaxis systems
	-m || --mfile              - output file with MCP chemoreceptros
	-b || --dbase              - specify database: mist or mist-mags
	'''

LOGGER = logging.getLogger(__name__)
logging.basicConfig(filename=sys.argv[0].replace(".py", "") + "_log.txt", level=logging.INFO)
TIMEOUT_FILE = sys.argv[0].replace(".py", "")  + "_timeout_info.txt"

INPUT_FILE = None
INPUT_FILE2 = None
OUTPUT_FILE_MST = None
OUTPUT_FILE_CHEM = None
OUTPUT_FILE_MCP = None

DATABASE = "mist"
GENOME_TO_PROTEIN_COUNT = {}

GENOMES_URL = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/"
METAGENOMES_URL = "https://metagenomes.asc.ohio-state.edu/v1/genomes/"
DATABASE_TO_URL = {"mist": GENOMES_URL, "mist-mags": METAGENOMES_URL}

MAJOR_ST_TO_COUNTS = OrderedDict([("total", 0), ("ocp_total", 0), ("tcp_total", 0), ("tcp_hk", 0), ("tcp_hhk", 0), ("tcp_rr", 0), ("tcp_hrr", 0), ("tcp_other", 0), ("ecf", 0), ("chem_sys", 0), ("other", 0), \
								 	("total_norm_by_protein_count", 0), ("ocp_total_norm_by_protein_count", 0), ("tcp_total_norm_by_protein_count", 0), ("tcp_hk_norm_by_protein_count", 0), \
									("tcp_hhk_norm_by_protein_count", 0), ("tcp_rr_norm_by_protein_count", 0), ("tcp_hrr_norm_by_protein_count", 0), ("tcp_other_norm_by_protein_count", 0), \
									("ecf_norm_by_protein_count", 0), ("chem_sys_norm_by_protein_count", 0), ("other_norm_by_protein_count", 0)])

CHEM_TEMPLATE_TO_COUNTS = OrderedDict([("protein_type", None), ("$total", 0), ("chew", 0), ("checx", 0), ("other", 0), ("MAC1", 0), ("MAC2", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("ACF", 0), ("Tfp", 0), \
								("$total_norm_by_protein_count", 0), ("chew_norm_by_protein_count", 0), ("checx_norm_by_protein_count", 0), ("other_norm_by_protein_count", 0), ("MAC1_norm_by_protein_count", 0), ("MAC2_norm_by_protein_count", 0), ("F1_norm_by_protein_count", 0), ("F2_norm_by_protein_count", 0), \
								("F3_norm_by_protein_count", 0), ("F4_norm_by_protein_count", 0), ("F5_norm_by_protein_count", 0), ("F6_norm_by_protein_count", 0), ("F7_norm_by_protein_count", 0), ("F8_norm_by_protein_count", 0), ("F9_norm_by_protein_count", 0), \
								("F10_norm_by_protein_count", 0), ("F11_norm_by_protein_count", 0), ("F12_norm_by_protein_count", 0), ("F13_norm_by_protein_count", 0), ("F14_norm_by_protein_count", 0), ("F15_norm_by_protein_count", 0), \
								("F16_norm_by_protein_count", 0), ("F17_norm_by_protein_count", 0), ("ACF_norm_by_protein_count", 0), ("Tfp_norm_by_protein_count", 0)])

CHEA_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEW_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHECX_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEMOTHER_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHED_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEZ_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEV_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()

CHEB_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHER_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
MCP_TO_COUNTS = OrderedDict([("$total", 0), ("64H", 0), ("58H", 0), ("52H", 0), ("48H", 0), ("44H", 0), ("42H", 0), ("40H", 0), ("38H", 0), ("36H", 0), ("34H", 0), ("28H", 0), ("24H", 0), \
								("total_norm_by_protein_count", 0), ("64H_norm_by_protein_count", 0), ("58H_norm_by_protein_count", 0), ("52H_norm_by_protein_count", 0), \
								("48H_norm_by_protein_count", 0), ("44H_norm_by_protein_count", 0), ("42H_norm_by_protein_count", 0), ("40H_norm_by_protein_count", 0), \
								("38H_norm_by_protein_count", 0), ("36H_norm_by_protein_count", 0), ("34H_norm_by_protein_count", 0), ("28H_norm_by_protein_count", 0), ("24H_norm_by_protein_count", 0)]) 

COMPONENT_TO_DATASTRUCTURE = dict([("chea", CHEA_TO_COUNTS), ("chew", CHEW_TO_COUNTS), ("checx", CHECX_TO_COUNTS), ("other", CHEMOTHER_TO_COUNTS), ("chev", CHEV_TO_COUNTS), ("cheb", CHEB_TO_COUNTS), ("cher", CHER_TO_COUNTS), ("ched", CHED_TO_COUNTS), ("chez", CHEZ_TO_COUNTS), ("mcp", MCP_TO_COUNTS)]) 

def initialize(argv):
	global INPUT_FILE, INPUT_FILE2, OUTPUT_FILE_MST, OUTPUT_FILE_MCP, OUTPUT_FILE_CHEM, DATABASE, GTDB_METADATA_FILE
	try:
		opts, args = getopt.getopt(argv[1:],"hi:s:o:c:m:b:g:",["help", "ifile=", "sfile=", "ofile=", "cfile=", "mfile=", "dbase=", "gtdbmeta="])
		if len(opts) == 0:
			raise getopt.GetoptError("Options are required\n")
	except getopt.GetoptError as e:
		print(("===========ERROR==========\n " + str(e) + USAGE))
		sys.exit(2)
	try:
		for opt, arg in opts:
			if opt in ("-h", "--help"):
				print(USAGE)
				sys.exit()
			elif opt in ("-i", "--ifile"):
				INPUT_FILE = str(arg).strip()
			elif opt in ("-s", "--sfile"):
				INPUT_FILE2 = str(arg).strip()
			elif opt in ("-o", "--ofile"):
				OUTPUT_FILE_MST = str(arg).strip()
			elif opt in ("-c", "--cfile"):
				OUTPUT_FILE_CHEM = str(arg).strip()
			elif opt in ("-m", "--mfile"):
				OUTPUT_FILE_MCP = str(arg).strip()
			elif opt in ("-b", "--dbase"):
				DATABASE = str(arg).strip()
				if DATABASE not in DATABASE_TO_URL:
					print ("Database should be one of the following: " + ", ".join(DATABASE_TO_URL.keys()))
					sys.exit(2)
			elif opt in ("-g", "--gtdbmeta"):
				GTDB_METADATA_FILE = str(arg).strip()
	except Exception as e:
		print(("===========ERROR==========\n " + str(e) + USAGE))
		sys.exit(2)
	
	with open(INPUT_FILE2, "r") as iFile2:
		for line in iFile2:
            #$0 - genome version, $1 - genome accession, $2 - genome size, $3 - protein counts, $4 - GTDB taxonomy, $5 - NCBI taxonomy
			record = line.strip().split("\t")
			genome_version = record[0]
			GENOME_TO_PROTEIN_COUNT[genome_version] = record[3]

def saveHeaders():
	# Save headers
	## Major ST systems
	with open(OUTPUT_FILE_MST, 'w') as ofile:
		ofile.write("genome\tgenome_accession\t" + "\t".join(list(MAJOR_ST_TO_COUNTS.keys())) + "\n")
	## Chemotaxis systems
	with open(OUTPUT_FILE_CHEM, 'w') as cfile, open(OUTPUT_FILE_MCP, 'w') as mfile:
		cfile.write("genome\tgenome_accession\t" + "\t".join(list(CHEM_TEMPLATE_TO_COUNTS.keys())).replace("$", "") + "\n")
		mfile.write("genome\tgenome_accession\t" + "\t".join(list(MCP_TO_COUNTS.keys())).replace("$", "") + "\n")

def collectSignalGenesCounts():
	recordCounter = 0
	with open(INPUT_FILE) as inputfile:
		for line in inputfile:
			genomeVersion = line.split("\t")[1]
			recordCounter+=1
			print((genomeVersion + "\t" + str(recordCounter)))
			constructedUrl = DATABASE_TO_URL[DATABASE] + genomeVersion + "/stp-matrix?per_page=1"
			for iteration in range(0, 10):
				try:
					result = urllib.request.urlopen(constructedUrl)
					data = json.loads(result.read().decode("utf-8"))
					if len(data) == 0:   #No data anymore from this page on
						break
					if "name" in data:	#404 NotFoundError
						break
				except (urllib.error.HTTPError, urllib.error.URLError, json.decoder.JSONDecodeError) as error:
					if iteration == 9:
						with open (TIMEOUT_FILE, "a") as timeoutFile:
							LOGGER.info("Ten attempts to retrieve data were unsuccessful. Save the genome caused the problem to %s file", TIMEOUT_FILE)
							timeoutFile.write(genomeVersion + "\n")
					#sleep 5 seconds if gateway timeout happened
					LOGGER.error("Timeout error: %s", error)
					LOGGER.info("Attempt " + str(iteration) + ". Sleep for 5 seconds...")
					time.sleep(5)
					LOGGER.info("Continue.")
					continue
		
				if "counts" in data:
					if "$total" in data["counts"]:
						MAJOR_ST_TO_COUNTS["total"] = data["counts"]["$total"]
					if "ocp" in data["counts"]:     
						MAJOR_ST_TO_COUNTS["ocp_total"] = data["counts"]["ocp"]["$total"]
					if "tcp" in data["counts"]: 
						MAJOR_ST_TO_COUNTS["tcp_total"] = data["counts"]["tcp"]["$total"]
						if "hk" in data["counts"]["tcp"]:
							MAJOR_ST_TO_COUNTS["tcp_hk"] = data["counts"]["tcp"]["hk"]["$total"]
						if "hhk" in data["counts"]["tcp"]:      
							MAJOR_ST_TO_COUNTS["tcp_hhk"] = data["counts"]["tcp"]["hhk"]["$total"]
						if "rr" in data["counts"]["tcp"]:       
							MAJOR_ST_TO_COUNTS["tcp_rr"] = data["counts"]["tcp"]["rr"]["$total"]
						if "hrr" in data["counts"]["tcp"]:      
							MAJOR_ST_TO_COUNTS["tcp_hrr"] = data["counts"]["tcp"]["hrr"]["$total"]
						if "other" in data["counts"]["tcp"]:            
							MAJOR_ST_TO_COUNTS["tcp_other"] = data["counts"]["tcp"]["other"]["$total"]
					if "chemotaxis" in data["counts"]:						
						if "chea" in data["counts"]["chemotaxis"]:
							#This is for the one file output with all ST proteins and the number of chemotaxis systems
							MAJOR_ST_TO_COUNTS["chem_sys"] = data["counts"]["chemotaxis"]["chea"]["$total"]						
						for component, data_strcuture in COMPONENT_TO_DATASTRUCTURE.items():
							processChemotaxisSystems(genomeVersion, data, component, data_strcuture)

					if "ecf" in data["counts"]:
						MAJOR_ST_TO_COUNTS["ecf"] = data["counts"]["ecf"]["$total"]
					if "other" in data["counts"]:
						MAJOR_ST_TO_COUNTS["other"] = data["counts"]["other"]["$total"]
				
				normalize_major_st_by_protein_counts(genomeVersion)

				# Save chemotaxis systems
				for component, dataDict in COMPONENT_TO_DATASTRUCTURE.items():
					if component != "mcp":
						saveToFile(dataDict, OUTPUT_FILE_CHEM, genomeVersion)
					else :
						saveToFile(dataDict, OUTPUT_FILE_MCP, genomeVersion)
					resetDataStructure(dataDict)

				# Save majsor ST systems
				saveToFile(MAJOR_ST_TO_COUNTS, OUTPUT_FILE_MST, genomeVersion)
				
				break

def processChemotaxisSystems(genomeVersion, data, entity, component_to_counts = False):
	postfix = "_norm_by_protein_count"
	component_to_counts["protein_type"] = entity
	if entity in data["counts"]["chemotaxis"]:
		for system, subdata in data["counts"]["chemotaxis"][entity].items():
			if system == "$total":
				component_to_counts[system] = subdata
				# normalization by protein count; the same below
				component_to_counts[system+postfix] = subdata/float(GENOME_TO_PROTEIN_COUNT[genomeVersion])
			elif system in component_to_counts:
				component_to_counts[system] = subdata["$total"]
				component_to_counts[system+postfix] = subdata["$total"]/float(GENOME_TO_PROTEIN_COUNT[genomeVersion])

def normalize_major_st_by_protein_counts(genomeVersion):
	systems = [system for system in MAJOR_ST_TO_COUNTS.keys() if not system.endswith("_norm_by_protein_count")]
	for sys in systems:
		MAJOR_ST_TO_COUNTS[sys+"_norm_by_protein_count"] = MAJOR_ST_TO_COUNTS[sys]/float(GENOME_TO_PROTEIN_COUNT[genomeVersion])

def  saveToFile(dataDict, outputFile, genomeVersion):
	with open(outputFile, 'a') as output_file:
		output_file.write(genomeVersion + "\t" + genomeVersion.split(".")[0] +  "\t" + "\t".join(map(str, (dataDict.values()))) + "\n")

def resetDataStructure(dataStrcture):
	for key in dataStrcture:
		dataStrcture[key] = 0

def main(argv):
	initialize(argv)
	saveHeaders()
	collectSignalGenesCounts()

main(sys.argv)

	








	

	
	
	
	
	
	
	
	
	
	
	
	
	
