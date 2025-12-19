#!/usr/bin/python3
import sys, getopt
import urllib.request, urllib.error
import json
import collections
import os.path
import time
import logging

USAGE = "\n\nThe script queries MiST db via it's API for signal transduction in genomes.\n\n" + \
	"It produces ST info for each genome. \n\n" + \
	"python 	" + sys.argv[0] + '''
	-h || --help               - help
	-i || --ifile              - input file (a new-line delimeted list of genome versions)
	-s || --sfile              - input file with ST-domain to classification information
	-o || --ofile              - output file
	-p || --pfile              - output file with genome id, protein id, and domain architecture
	-t || --task               - what to do. One of the following:
								 identifyDomainCombinsAcrossComps, identifyDomainCombinsAcrossGenomes; classifySignalGenesAcrossCompons, classifySignalGenesInCompons
	-b || --dbase              - specify database: mist or mist-mags
	'''

LOGGER = logging.getLogger(__name__)
logging.basicConfig(filename=sys.argv[0].replace(".py", "") + "_log.txt", level=logging.INFO)
TIMEOUT_FILE = sys.argv[0].replace(".py", "")  + "_timeout_info.txt"

FUNCTIONAL_CATEGORIES = collections.OrderedDict([("TR", 0), ("K/P", 0), ("PP",0), ("UNK", 0), ("DNA_UNK", 0), ("Che", 0)])


HEADER_TO_INDEX_IN_ST_FILE = \
	{"Classification_marker": 0, "Definition":1, "Domain_name":2, "Pfam_clan_(superfamily)":3, \
	"Description":4, "Function":5, "Specific_to_signal_transduction":6, "Source":7, "New_classification":8, \
	"Abbr_new_classification":9}

GENOME_VERSIONS = list()
DOMAIN_NAME_TO_NEW_CLASSIFICATION = dict()


# GENOMES_URL = "https://api.mistdb.caltech.edu/v1/genomes/"
# #COMPONENT_INFO_FIELDS = "?fields=id,version&fields.Component=id,name,version,length,type"
# SIGNAL_GENES_ADDITIONAL_FIELDS = "/signal-genes?count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"

DATABASE = "mist"

GENOMES_URL = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/"
METAGENOMES_URL = "https://metagenomes.asc.ohio-state.edu/v1/genomes/"
DATABASE_TO_URL = {"mist": GENOMES_URL, "mist-mags": METAGENOMES_URL}
SIGNAL_GENES_ADDITIONAL_FIELDS = "/signal-genes?count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"

INPUT_FILE = None
ST_INPUT_FILE = None
OUTPUT_FILE = None
OUTPUT_FILE2 = None
TASK = None

def initialize(argv):
	global INPUT_FILE, OUTPUT_FILE, OUTPUT_FILE2, ST_INPUT_FILE, TASK, DATABASE
	try:
		opts, args = getopt.getopt(argv[1:],"hi:s:o:p:t:b:",["help", "ifile=", "sfile=", "ofile=", "pfile=", "task=", "dbase="])
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
			elif opt in ("-s", "--sfile"):
				ST_INPUT_FILE = str(arg).strip()
			elif opt in ("-o", "--ofile"):
				OUTPUT_FILE = str(arg).strip()
			elif opt in ("-p", "--pfile"):
				OUTPUT_FILE2 = str(arg).strip()
			elif opt in ("-t", "--task"):
				TASK = str(arg).strip()
			elif opt in ("-b", "--dbase"):
				DATABASE = str(arg).strip()
				if DATABASE not in DATABASE_TO_URL:
					print ("Database should be one of the following: " + ", ".join(DATABASE_TO_URL.keys()))
					sys.exit(2)
	except Exception as e:
		print("===========ERROR==========\n " + str(e) + USAGE)
		sys.exit(2)

def initialyzeSTCollectionAndGenomeList():
	with open(ST_INPUT_FILE, "r") as stInputFile:
		for record in stInputFile:
			record = record.split("\t")
			domainColumn = HEADER_TO_INDEX_IN_ST_FILE["Domain_name"]
			classificColumn = HEADER_TO_INDEX_IN_ST_FILE["Abbr_new_classification"]
			classification = record[classificColumn].strip()
			if len(classification):
				DOMAIN_NAME_TO_NEW_CLASSIFICATION[record[domainColumn].strip()] = classification
	global GENOME_VERSIONS
	with open(INPUT_FILE, "r") as inputFile:
		GENOME_VERSIONS = [record.strip() for record in inputFile]

def retrieveSignalGenesFromMist(genomeVersion):
	#genomeURL = GENOMES_URL + genomeVersion
	genomeURL = DATABASE_TO_URL[DATABASE] + genomeVersion
	signalGenesList = list()
	noDataAnymore = False
	for num in range(1, 101):
		for iteration in range(0, 10):
			try:
				additionaFields = SIGNAL_GENES_ADDITIONAL_FIELDS.replace("%PAGE%", str(num))
				constructedUrl = genomeURL + additionaFields 
				result = urllib.request.urlopen(constructedUrl)
				resultAsJson = json.loads(result.read().decode("utf-8"))
				if len(resultAsJson) == 0:   #No data anymore from this page on
					noDataAnymore = True
					break
				if "name" in resultAsJson:	#404 NotFoundError
					break
			# #except ValueError: #504 Gateway timeouts
			# except urllib.error.HTTPError:
			# #except json.decoder.JSONDecodeError:   #504 Gateway timeouts  From Python 3.5+
			except (urllib.error.HTTPError, urllib.error.URLError, json.decoder.JSONDecodeError) as error:
				if iteration == 9:
					# with open (TIMEOUT_FILE, "a") as timeoutFile:
					# 	timeoutFile.write(genomeVersion + "\n")
					with open (TIMEOUT_FILE, "a") as timeoutFile:
						LOGGER.info("Ten attempts to retrieve data were unsuccessful. Save the genome caused the problem to %s file", TIMEOUT_FILE)
						timeoutFile.write(genomeVersion + "\n")
				#sleep 5 seconds if gateway timeout happened
				LOGGER.error("Timeout error: %s", error)
				LOGGER.info("Attempt " + str(iteration) + ". Sleep for 5 seconds...")
				time.sleep(5)
				LOGGER.info("Continue.")
				continue
			signalGenesList.extend(resultAsJson)
			break
		if noDataAnymore:
			break
	return signalGenesList			
	
def processSignalGenes():
	counter = 1
	if TASK == "classifySignalGenesInCompons":
		if not os.path.exists(OUTPUT_FILE):
			with open(OUTPUT_FILE, "w") as outputFile:
				outputFile.write("Genome_Version\tComponent_Name\tComponent_Version\t" + "\t".join(list(FUNCTIONAL_CATEGORIES.keys())) + "\n")
	elif TASK == "classifySignalGenesAcrossCompons":
		if not os.path.exists(OUTPUT_FILE):
			with open(OUTPUT_FILE, "w") as outputFile:
				outputFile.write("Genome_Version\t" + "\t".join(list(FUNCTIONAL_CATEGORIES.keys())) + "\n")
	for genomeVersion in GENOME_VERSIONS:
		if genomeVersion.strip():
			print(genomeVersion + "\t" + str(counter))
			counter+=1
			signalGenesList = retrieveSignalGenesFromMist(genomeVersion)
			if TASK == "classifySignalGenesInCompons":
				classifySignalGenesInCompons(genomeVersion, signalGenesList)
			elif TASK == "classifySignalGenesAcrossCompons":
				classifySignalGenesAcrossCompons(genomeVersion, signalGenesList)
		
def classifySignalGenesInCompons(genomeVersion, signalGenesList,):
	componentToFunctionalCategories = dict()
	for gene in signalGenesList:
		componentNameAndVersion = gene["Component"]["name"] + "\t" + gene["Component"]["version"]
		if "counts" in gene:
			for domain in gene["counts"]:
				domain = domain.strip()
				if domain in DOMAIN_NAME_TO_NEW_CLASSIFICATION:
					if componentNameAndVersion not in componentToFunctionalCategories:
						componentToFunctionalCategories[componentNameAndVersion] = FUNCTIONAL_CATEGORIES.copy() 
					componentToFunctionalCategories[componentNameAndVersion][DOMAIN_NAME_TO_NEW_CLASSIFICATION[domain]]+=1		
	
	for componentNameAndVersion, categories in list(componentToFunctionalCategories.items()):
		resultString = genomeVersion + "\t" + componentNameAndVersion
		for category, count in list(categories.items()):
			resultString = resultString + "\t" + str(count)
		with open(OUTPUT_FILE, "a") as outputFile:
			outputFile.write(resultString + "\n")

def classifySignalGenesAcrossCompons(genomeVersion, signal_genes_list):
	functional_category_to_count = FUNCTIONAL_CATEGORIES.copy()
	resultString = genomeVersion
	for gene in signal_genes_list:
		if "counts" in gene:
			for domain in gene["counts"]:
				domain = domain.strip()
				if domain in DOMAIN_NAME_TO_NEW_CLASSIFICATION:
					functional_category_to_count[DOMAIN_NAME_TO_NEW_CLASSIFICATION[domain]]+=1
				
	for category, count in list(functional_category_to_count.items()):
		resultString = resultString + "\t" + str(count)
	
	with open(OUTPUT_FILE, "a") as outputFile:
		outputFile.write(resultString + "\n")
	

##*********************************************************************##
##********************** Domains processing block**********************##
def processDomains():
	if TASK == "identifyDomainCombinsAcrossComps":
		genomeNumber = 1
		for genomeVersion in GENOME_VERSIONS:
			if genomeVersion.strip():
				print("Genome Number: " + str(genomeNumber))
				genomeNumber+=1
				signalGenesList = retrieveSignalGenesFromMist(genomeVersion)
				domainCombinToCount = collections.defaultdict(int)
				for gene in signalGenesList:
					dominsListSordtedAsStr = prepareDomains(gene, domainCombinToCount)
					if dominsListSordtedAsStr:
						with open(OUTPUT_FILE2, "a") as outputFile:
							outputFile.write(genomeVersion + "\t" + gene["Gene"]["version"] + "\t" + dominsListSordtedAsStr + "\n")
						
				domainCombinToCountList = sorted(list(domainCombinToCount.items()), key=lambda a: a[1], reverse=True)
				with open(OUTPUT_FILE, "a") as outputFile:
					for domainCombinAndCount in domainCombinToCountList:
						outputFile.write(genomeVersion + "\t" + domainCombinAndCount[0] + "\t" + str(domainCombinAndCount[1]) + "\n")	

						
	elif TASK == "identifyDomainCombinsAcrossGenomes":
		domainCombinAcrossGenomesToCount = collections.defaultdict(int) 
		for genomeVersion in GENOME_VERSIONS:
			if genomeVersion.strip():
				signalGenesList = retrieveSignalGenesFromMist(genomeVersion) 
				for gene in signalGenesList:
					prepareDomains(gene, domainCombinAcrossGenomesToCount)
		domainCombinAndCountList = sorted(list(domainCombinAcrossGenomesToCount.items()), key=lambda a: a[1], reverse=True)
		with open(OUTPUT_FILE, "w") as outputFile:
			for domainCombinAndCount in domainCombinAndCountList:
				outputFile.write(domainCombinAndCount[0] + "\t" + str(domainCombinAndCount[1]) + "\n")			

def prepareDomains(gene, domainCombinToCount):
	if "Gene" in gene and "Aseq" in gene["Gene"] and "pfam31" in gene["Gene"]["Aseq"]:
		domainsSorted = sorted(gene["Gene"]["Aseq"]["pfam31"], key=lambda x: x["ali_from"], reverse=False)
		if len (domainsSorted) > 0:
			domainsFiltered = removeOverlapps(domainsSorted)
			#domainsFilteredNames = [domain["name"] for domain in domainsFiltered]
			domainsFilteredNames = set()
			for domain in domainsFiltered:
				domainFinal = domain["name"].strip()
				domainSplitted = domain["name"].split("_")
				if domainSplitted[0] == "PAS" or domainSplitted[0] == "Sigma70":
					domainFinal = domainSplitted[0]
				domainsFilteredNames.add(domainFinal)	
			dominsListSordtedAsStr = ",".join(sorted(domainsFilteredNames))
			domainCombinToCount[dominsListSordtedAsStr]+=1
			return dominsListSordtedAsStr
	return False
					
		
def removeOverlapps(domainsSorted):
	tolerance = 10
	pfam31Final = []
	pfam1 = domainsSorted[0]
	significantPfam = pfam1
	overlapLength = None
	lastAdded = pfam1
	pfam31Final.append(pfam1)
	for pfam2 in domainsSorted[1:]:
		if pfam1["ali_to"] > pfam2["ali_from"]:
			overlapLength = pfam1["ali_to"] - pfam2["ali_from"]
			if overlapLength > tolerance:
				significantPfam = compareEvalues(pfam1, pfam2)
				# if the previously added is pfam1 and it's less significant than pfam 2
				# then remove this previously added and add pfam 2
				if lastAdded == pfam1 and lastAdded != significantPfam:
					pfam31Final.remove(lastAdded)
					pfam31Final.append(significantPfam)
				lastAdded = significantPfam
			else:
				pfam31Final.append(pfam2)
				lastAdded = pfam2
				significantPfam = pfam2
			pfam1 = significantPfam
		else:
			pfam31Final.append(pfam2)
			lastAdded = pfam2
			significantPfam = pfam2
			pfam1 = significantPfam
	return pfam31Final

def compareEvalues(pfam1, pfam2):
	if "i_evalue" in pfam1:
		eval1 = pfam1["i_evalue"]
		eval2 = pfam2["i_evalue"]
	significantPfam = None
	if eval1 < eval2:
		significantPfam = pfam1
	elif eval1 > eval2:
		significantPfam = pfam2
	elif eval1 == eval2:
		if (pfam1["ali_to"] - pfam1["ali_from"]) >= (pfam2["ali_to"] - pfam2["ali_from"]):
			significantPfam = pfam1
		else:
			significantPfam = pfam2
	return significantPfam

##***************** Domains processing block finish *******************##		
##*********************************************************************##		
		
def main(argv):
	initialize(argv)
	initialyzeSTCollectionAndGenomeList()
	if TASK == "classifySignalGenesAcrossCompons" or TASK == "classifySignalGenesInCompons":
		print ("Identifying and processing Signal genes")
		processSignalGenes()
	elif TASK == "identifyDomainCombinsAcrossComps" or TASK == "identifyDomainCombinsAcrossGenomes":
		print ("Identifying and processing Signal domains")
		processDomains()

	

main(sys.argv)
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
