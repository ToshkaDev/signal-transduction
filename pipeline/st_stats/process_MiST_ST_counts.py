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
	-o || --ofile              - output file
	-b || --dbase              - specify database: mist or mist-mags
	-g || --gtdbmeta           - GTDB metadata filed prepared by the SIGNAL pipeline 
	'''

LOGGER = logging.getLogger(__name__)
logging.basicConfig(filename=sys.argv[0].replace(".py", "") + "_log.txt", level=logging.INFO)
TIMEOUT_FILE = sys.argv[0].replace(".py", "")  + "_timeout_info.txt"

INPUT_FILE = None
OUTPUT_FILE = None

DATABASE = "mist"
GTDB_METADATA_FILE = "../../input/gtdb_metadata/ar_bac_metadata_r214_p.tsv"

GENOMES_URL = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/"
METAGENOMES_URL = "https://metagenomes.asc.ohio-state.edu/v1/genomes/"
DATABASE_TO_URL = {"mist": GENOMES_URL, "mist-mags": METAGENOMES_URL}

CHEA_FILE_OUTPUT = sys.argv[0].replace(".py", "") + "_CheA_across_genomes.txt"
CHEB_FILE_OUTPUT = sys.argv[0].replace(".py", "") + "_CheB_across_genomes.txt"
CHER_FILE_OUTPUT = sys.argv[0].replace(".py", "") + "_CheR_across_genomes.txt"
CHEZ_FILE_OUTPUT = sys.argv[0].replace(".py", "") + "_CheZ_across_genomes.txt"
CHED_FILE_OUTPUT = sys.argv[0].replace(".py", "") + "_CheD_across_genomes.txt"
CHEV_FILE_OUTPUT = sys.argv[0].replace(".py", "") + "_CheV_across_genomes.txt"
MCP_FILE_OUTPUT = sys.argv[0].replace(".py", "") + "_MCP_across_genomes.txt"

TEMPLATE_TO_COUNTS = OrderedDict([("$total", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0)])
CHEA_TO_COUNTS = OrderedDict([("$total", 0), ("chew", 0), ("checx", 0), ("other", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0)])
CHED_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEZ_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEV_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEB_TO_COUNTS = OrderedDict([("$total", 0), ("MAC1", 0), ("MAC2", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0)])
CHER_TO_COUNTS = CHEB_TO_COUNTS.copy()
MCP_TO_COUNTS = OrderedDict([("$total", 0), ("64H", 0), ("58H", 0), ("52H", 0), ("48H", 0), ("44H", 0), ("42H", 0), ("40H", 0), ("38H", 0), ("36H", 0), ("34H", 0), ("28H", 0), ("24H", 0)]) 

OUTPU_FILE_TO_DICT = OrderedDict([(CHEA_FILE_OUTPUT, CHEA_TO_COUNTS), (CHEB_FILE_OUTPUT, CHEB_TO_COUNTS), (CHER_FILE_OUTPUT, CHER_TO_COUNTS), (CHEZ_FILE_OUTPUT, CHEZ_TO_COUNTS), (CHED_FILE_OUTPUT, CHED_TO_COUNTS), (CHEV_FILE_OUTPUT, CHEV_TO_COUNTS), (MCP_FILE_OUTPUT, MCP_TO_COUNTS) ])
COMPONENT_TO_DATASTRUCTURE = dict([("chea", CHEA_TO_COUNTS), ("chew", CHEA_TO_COUNTS), ("checx", CHEA_TO_COUNTS), ("other", CHEA_TO_COUNTS), ("chev", CHEV_TO_COUNTS), ("cheb", CHEB_TO_COUNTS), ("cher", CHER_TO_COUNTS), ("ched", CHED_TO_COUNTS), ("chez", CHEZ_TO_COUNTS), ("mcp", MCP_TO_COUNTS)]) 

GENOME_TO_TAXONOMY = {}

def initialize(argv):
	global INPUT_FILE, OUTPUT_FILE, TASK, DATABASE, GTDB_METADATA_FILE
	try:
		opts, args = getopt.getopt(argv[1:],"hi:o:b:g:",["help", "ifile=", "ofile=", "dbase=", "gtdbmeta"])
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
			elif opt in ("-o", "--ofile"):
				OUTPUT_FILE = str(arg).strip()
				TASK = str(arg).strip()
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
	
	with open(GTDB_METADATA_FILE) as gtdbMetaFile:
		for line in gtdbMetaFile:
			line_s = line.split("\t")
			GENOME_TO_TAXONOMY[line_s[0]] = line_s[4]


def collectSignalGenesCounts():
	with open(OUTPUT_FILE, 'w') as output_file:
		output_file.write( "\t".join('genomeID, Total, OCP_total, TCP_total, TCP_HK, TCP_HHK, TCP_RR, TCP_HRR, TCP_Other, ChemSys, ECF, Other, Taxonomy, \n'.split(",")) )
	for outpufFile, data_structure in OUTPU_FILE_TO_DICT.items():
		with open(outpufFile, 'w') as outpuf_file:
			outpuf_file.write("genomeID\t" + "\t".join(data_structure.keys()) + "\t" + "Taxonomy\n")
	
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
		
				total = OCP_total = TCP_total = TCP_HK = TCP_HHK = TCP_RR = TCP_HRR = TCP_Other = ECF = ChemSys = Other = "0"
				if "counts" in data:
					if "$total" in data["counts"]:
						total = data["counts"]["$total"]
					if "ocp" in data["counts"]:     
						OCP_total = data["counts"]["ocp"]["$total"]
					if "tcp" in data["counts"]: 
						TCP_total = data["counts"]["tcp"]["$total"]
						if "hk" in data["counts"]["tcp"]:
							TCP_HK = data["counts"]["tcp"]["hk"]["$total"]
						if "hhk" in data["counts"]["tcp"]:      
							TCP_HHK = data["counts"]["tcp"]["hhk"]["$total"]
						if "rr" in data["counts"]["tcp"]:       
							TCP_RR = data["counts"]["tcp"]["rr"]["$total"]
						if "hrr" in data["counts"]["tcp"]:      
							TCP_HRR = data["counts"]["tcp"]["hrr"]["$total"]
						if "other" in data["counts"]["tcp"]:            
							TCP_Other = data["counts"]["tcp"]["other"]["$total"]
					if "chemotaxis" in data["counts"]:						
						if "chea" in data["counts"]["chemotaxis"]:
							#This is for the one file output with all ST proteins and the number of chemotaxis systems
							ChemSys = data["counts"]["chemotaxis"]["chea"]["$total"]						
						for component, data_strcuture in COMPONENT_TO_DATASTRUCTURE.items():
							processChemotaxisSystems(data, component, data_strcuture)

					if "ecf" in data["counts"]:
						ECF = data["counts"]["ecf"]["$total"]
					if "other" in data["counts"]:
						Other = data["counts"]["other"]["$total"]
				
				for outputFile, dataDict in OUTPU_FILE_TO_DICT.items():
					saveToFile(dataDict, outputFile, genomeVersion)
					resetDataStructure(dataDict)
				with open(OUTPUT_FILE, 'a') as output_file:
					output_file.write( "\t".join(map(str, ([genomeVersion, total, OCP_total, TCP_total, TCP_HK, TCP_HHK, TCP_RR, TCP_HRR, TCP_Other, ChemSys, ECF, Other, GENOME_TO_TAXONOMY[genomeVersion]]))) + "\n") 
				 	
				break

def processChemotaxisSystems(data, entity, component_to_counts = False):
	#entity can be chea, mcp, cheb, cher, and so on
	additional = set(["checx", "chew", "other"])
	if entity in data["counts"]["chemotaxis"]:
		for system, subdata in data["counts"]["chemotaxis"][entity].items():
			if system == "$total":
				if entity not in additional:
					component_to_counts[system] = subdata
				else:
					component_to_counts[entity] = subdata
			elif system in component_to_counts:
				component_to_counts[system] = subdata["$total"]	

def  saveToFile(dataDict, outputFile, genomeVersion):
	with open(outputFile, 'a') as output_file:
		output_file.write(genomeVersion + "\t" + "\t".join(map(str, (dataDict.values()))) + "\t" +  GENOME_TO_TAXONOMY[genomeVersion] + "\n")

def resetDataStructure(dataStrcture):
	for key in dataStrcture:
		dataStrcture[key] = 0

def main(argv):
	initialize(argv)
	collectSignalGenesCounts()

	
main(sys.argv)

	








	

	
	
	
	
	
	
	
	
	
	
	
	
	
