#!/usr/bin/python3
import sys, getopt
import urllib.request, urllib.parse, urllib.error
import json
import collections
import os.path
import time
import logging

OUT_FILE_HEADERS = ["Genome_id", "NCBI_id", "MiST_id", "protein_length", "domain_architecture", "sensors_or_regulators", "domain_counts", "domain_combinations", "\n"]

USAGE = "\n\nThe script queries MiST db via it's API for histidine kinases and response regulators in genomes. \n" + \
	"It outputs complete domain information and other data in tabulated format for both histidine kinases and response regualtors separately. \n" + \
	"Output fields: " + ", ".join(OUT_FILE_HEADERS).rstrip(", \n") + " \n" + \
	"python 	" + sys.argv[0] + '''
	-h || --help               - help
	-i || --ifile              - input file
	-f || --ffile              - first output file
	-s || --sfile              - second output file
	-d || --database           - specify database: mist or mist-mags
	-c || --continue           - start a new analysis or continue with allready existing provided files.
	                             Users are simply expected to specify -c (--continue) without provinding arguments.
	                             Default is without this paraeter specified, i.e. start a new analysis.
'''

#Variables controlled by the script parameters
INPUT_FILE = None
OUTPUT_FILE1 = "output_HK.tsv"
OUTPUT_FILE2 = "output_RR.tsv"
CONTINUE = False

#Variables set within the script
PROTEIN_TYPES = ["sensKinase", "respReg"]
PROTEIN_TYPE_TO_OUTFILE = {PROTEIN_TYPES[0]: OUTPUT_FILE1, PROTEIN_TYPES[1]: OUTPUT_FILE2}
GENOME_VERSIONS = None
TIMEOUT_FILE = "timeout_genomes.txt"
DATABASE = "mist"
LOGGER = logging.getLogger(__name__)
logging.basicConfig(filename=sys.argv[0].replace(".py", "") + "_log.txt", level=logging.INFO)

GENOMES_URL = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/"
METAGENOMES_URL = "https://metagenomes.asc.ohio-state.edu/v1/genomes/"
DATABASE_TO_URL = {"mist": GENOMES_URL, "mist-mags": METAGENOMES_URL}

STP_MATRIX = "/stp-matrix?page=%PAGE%&per_page=100"
SIGNAL_GENES_HK = "/signal-genes?where.component_id=%COMPONENT_ID%&where.ranks=tcp,hk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"
SIGNAL_GENES_HHK = "/signal-genes?where.component_id=%COMPONENT_ID%&where.ranks=tcp,hhk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"
SIGNAL_GENES_RR = "/signal-genes?where.component_id=%COMPONENT_ID%&where.ranks=tcp,rr&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"
SIGNAL_GENES_HRR = "/signal-genes?where.component_id=%COMPONENT_ID%&where.ranks=tcp,hrr&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"

HIS_KINASE_DIM_DOMAINS = ["HisKA", "HisKA_2", "HisKA_3", "H-kinase_dim", "His_kinase"]
HIS_KINASE_CATAL_DOMAINS = ["HATPase_c", "HATPase_c_2", "HATPase_c_5", "HWE_HK"]
RESPONSE_REG_DOMAINS = ["Response_reg", "FleQ"]

def initialize(argv):
	global INPUT_FILE, OUTPUT_FILE1, OUTPUT_FILE2, GENOME_VERSIONS, PROTEIN_TYPE_TO_OUTFILE, DATABASE, CONTINUE
	try:
		opts, args = getopt.getopt(argv[1:],"hi:f:s:d:c",["help", "ifile=", "ffile=", "sfile=", "database=", "continue"])
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
			elif opt in ("-f", "--ffile"):
				OUTPUT_FILE1 = str(arg).strip()
			elif opt in ("-s", "--sfile"):
				OUTPUT_FILE2 = str(arg).strip()
			elif opt in ("-d", "--database"):
				DATABASE = str(arg).strip()
				if DATABASE not in DATABASE_TO_URL:
					print ("Database should be one of the following: " + ", ".join(DATABASE_TO_URL.keys()))
					sys.exit(2)
			elif opt in ("-c", "--continue"):
				CONTINUE = True
	except Exception as e:
		print("===========ERROR==========\n " + str(e) + USAGE)
		sys.exit(2)
	#Initialize the dictionary with the provided files
	PROTEIN_TYPE_TO_OUTFILE = {PROTEIN_TYPES[0]: OUTPUT_FILE1, PROTEIN_TYPES[1]: OUTPUT_FILE2}
	if not CONTINUE:
		#Create ouput files and write headers:
		for oFile in PROTEIN_TYPE_TO_OUTFILE.values():
			with open(oFile, "w") as outFile:
				outFile.write("\t".join(OUT_FILE_HEADERS))

#The function first retreives stp-matrix and after analyzing the matrix it retrives signal genes for those components that have target signaling systems
#This is done, because it is much faster this way
def retrieveSignalGenesFromMist(genomeVersion, additionaFieldsTemplate):
	sensorRegulatorType = "hk"
	if additionaFieldsTemplate == SIGNAL_GENES_HK:
		sensorRegulatorType = "hk"
	elif additionaFieldsTemplate == SIGNAL_GENES_HHK:
		sensorRegulatorType = "hhk"
	elif additionaFieldsTemplate == SIGNAL_GENES_RR:
		sensorRegulatorType = "rr"
	elif additionaFieldsTemplate == SIGNAL_GENES_HRR:
		sensorRegulatorType = "hrr"
	genomeURL = DATABASE_TO_URL[DATABASE] + genomeVersion

	#Get stp-matrix and look at the numnber of components and save those components that have two-component systems.
	#They will be saved in componentsWithTcp list.
	componentsWithTcp = list() #componentsWithTcp will be populated
	getSignalGenes(genomeURL + STP_MATRIX, componentsWithTcp, genomeVersion, sensorRegulatorType, False, False)

	#Retrieve signal genes in those genomic components (chromosomes, scaffolds, or contigs depending on the assembly level) that have two-component systems
	signalGeneList = list()
	for component in componentsWithTcp:
		getSignalGenes(genomeURL, signalGeneList, genomeVersion, False, additionaFieldsTemplate, component)

	return signalGeneList			

def getSignalGenes(url, elementList, genomeVersion, tcpMatrix, additionaFieldsTemplate=False, component=False):
	noDataAnymore = False
	#This range is enough as 100 requests for 100 items per page is 10000 items. The numer of components per genome (or histidine kinase|respose regulators)
	#is an order of magnitde less than this.
	#'break' statement will stop sending requests when no data could be retrieved anymore
	for num in range(1, 101):
		if additionaFieldsTemplate:
			additionaFields = additionaFieldsTemplate.replace("%COMPONENT_ID%", str(component["id"])).replace("%PAGE%", str(num))
			constructedUrl = url + additionaFields
		else:
			constructedUrl = url.replace("%PAGE%", str(num))
		noDataAnymore = signalGenesRetriever(constructedUrl, elementList, genomeVersion, tcpMatrix, noDataAnymore)
		if noDataAnymore:
			break

def signalGenesRetriever(url, elementList, genomeVersion, tcpMatrix, noDataAnymore):
	for iteration in range (1, 11):
		try:
			result = urllib.request.urlopen(url)
			resultAsJson = json.loads(result.read().decode("utf-8"))
			#In case of tcpMatrix: No data anymore from this page on
			if tcpMatrix and "components" in resultAsJson and not resultAsJson["components"]:
				noDataAnymore = True
				break
			#Regular case of retrieveing proteins
			#No data anymore from this page on
			elif not resultAsJson:
				noDataAnymore = True
				break
			#404 NotFoundError
			if "name" in resultAsJson:
				break
		except (urllib.error.HTTPError, urllib.error.URLError) as error:
			if iteration == 10:
				with open (TIMEOUT_FILE, "a") as timeoutFile:
					LOGGER.info("Ten attempts to retrieve data were unsuccessful. Save the genome caused the problem to %s file", TIMEOUT_FILE)
					timeoutFile.write(genomeVersion + "\n")
			#sleep 5 seconds if gateway timeout happened
			LOGGER.error("Timeout error: %s", error)
			LOGGER.info("Attempt " + str(iteration) + ". Sleep for 5 seconds...")
			time.sleep(5)
			LOGGER.info("Continue.")
			continue
			
		if tcpMatrix:
			if "tcp" in resultAsJson["counts"]:
				for component in resultAsJson["components"]:
					if "tcp" in component["counts"] and tcpMatrix in component["counts"]["tcp"]:
						elementList.append(component)
		else:
			elementList.extend(resultAsJson)
		break
	return noDataAnymore

def processDomains():
	genomeNumber = 1
	with open(INPUT_FILE, "r") as inputFile:
		for genomeVersion in inputFile:
			genomeVersion = genomeVersion.split("\t")[1]
			if genomeVersion:
				print(" ".join(["Genome Number:", str(genomeNumber), "   Genome ID:", genomeVersion]))
				genomeNumber+=1
				listOfSignalGeneLists = []
				listOfSignalGeneLists.append((retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_HK), PROTEIN_TYPES[0]))
				listOfSignalGeneLists.append((retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_HHK), PROTEIN_TYPES[0]))
				listOfSignalGeneLists.append((retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_RR), PROTEIN_TYPES[1]))
				listOfSignalGeneLists.append((retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_HRR), PROTEIN_TYPES[1]))
				for signalGeneList in listOfSignalGeneLists:
					for gene in signalGeneList[0]:
						prepareDomains(gene, genomeVersion, signalGeneList[1])

##*********************************************************************##
##********************** Domains processing block**********************##
def prepareDomains(gene, genomeVersion, proteinType):
	if "Gene" in gene and "Aseq" in gene["Gene"] and "pfam31" in gene["Gene"]["Aseq"]:
		#Ordering domains according to how they are encoded in the gene
		domainsSorted = sorted(gene["Gene"]["Aseq"]["pfam31"], key=lambda x: x["ali_from"], reverse=False)
		#if there are entities in domainsSorted
		#if domainsSorted:
		if len (domainsSorted) > 0:
			domainsFiltered = removeOverlapps(domainsSorted)
			
			refseqVersion = gene["Gene"]["version"]
			geneStableId = gene["Gene"]["stable_id"]
			proteinLength = str(int(gene["Gene"]["length"]/3) - 1)
			domainsOutput = processHoles(domainsFiltered)

			domainArchitecture = ""
			domainArchitectureSensOrRegDomsOnly = []
			for domain in domainsOutput:
				domainArchitecture = domainArchitecture + "{}:{}-{},".format(domain["name"], domain["env_from"], domain["env_to"])
				domainName = domain["name"].lstrip("<").rstrip(">")
				if domain["name"] != "hole" and domainName not in HIS_KINASE_CATAL_DOMAINS and domainName not in HIS_KINASE_DIM_DOMAINS and domainName not in RESPONSE_REG_DOMAINS:
					domainArchitectureSensOrRegDomsOnly.append(domain["name"])

			#Generate a set of unique domain names and domain to count uniformly sorted
			#{'domain1': 1, 'domain2': 2}
			domainToCount = collections.defaultdict(int)
			for domain in domainsOutput:
				domainToCount[domain["name"]]+=1
			sortedDomainNames = sorted(domainToCount.keys())
			domainsFilteredNamesUniqueStr = ",".join(sortedDomainNames)
			domainsFilteredNamesUniqueCountsStr = ",".join(["{}:{}".format(domain, domainToCount[domain]) for domain in sortedDomainNames])

			#Save everything to a file:
			with open(PROTEIN_TYPE_TO_OUTFILE[proteinType], "a") as outputFile:
				outputRecord = "\t".join([genomeVersion, refseqVersion, geneStableId, proteinLength, domainArchitecture.rstrip(","), ",".join(domainArchitectureSensOrRegDomsOnly), domainsFilteredNamesUniqueCountsStr, domainsFilteredNamesUniqueStr])
				outputFile.write(outputRecord + "\n")

#### Process holes and domains BEGIN ####
def processHoles(domains):
	minDomainLength = 100
	minLengthForHisKA = 150
	domainsOutput = []
	HisKA_ind = -1
	HATPase_ind = -1
	for domain in domains:
		if domain["name"] in HIS_KINASE_DIM_DOMAINS: HisKA_ind = domains.index(domain)
		elif domain["name"] in HIS_KINASE_CATAL_DOMAINS: HATPase_ind = domains.index(domain)

	if HATPase_ind >= 0:
		if domains[:HATPase_ind]:
			if HisKA_ind >=0:
				checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength)
			#if HisKA domain has not been rcognized
			else:
				checkAndAddHolesAndDomains(domains[:HATPase_ind], domainsOutput, minDomainLength)
				#process the hole between the penultimate domain and the HATPase domain under the assumption that between these two domains the HisKA domain can be present
				HisKAprocessing(domains[HATPase_ind], domains[HATPase_ind-1], minLengthForHisKA, minDomainLength, domainsOutput)
				if domains[HATPase_ind+1:]:
					checkAndAddHolesAndDomains(domains[HATPase_ind+1:], domainsOutput, minDomainLength, domains[HATPase_ind]["env_to"])
		else:
			#if no domains upstream of the HATPase domain
			checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength)
	#if HisKA is present but HATPase is not (as HATPase cases are processed above)
	elif HisKA_ind >=0:
		checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength)
	#if both HisKA and HATPase are not in my two lists, HIS_KINASE_DIM_DOMAINS and HIS_KINASE_CATAL_DOMAINS, still process all the domains:
	#or if the domains belonw to response regulator, do this:
	else:
		checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength)
		
	return domainsOutput

def HisKAprocessing(domainUltimate, domainPenultimate, minLengthForHisKA, minDomainLength, domainsOutput):
	#process the hole between the penultimate domain and the HATPase domain, under the assumption that between these two domains the HisKA domain can be present
	#            domainPenultimate domainUltimate
	#domains: d1, <hole>, <HisKA>, HTPAse_c
	HisKAspace = domainUltimate["env_from"] - domainPenultimate["env_to"]        
	holeEnd = domainPenultimate["env_to"] + (HisKAspace - minLengthForHisKA)
	HisKAdomain = {"name": "<HisKA>", "env_from": holeEnd+1, "env_to": domainUltimate["env_from"]-1}                          																						#domains[-2]         domains[-1] 
	#if the below condition is satisfied then there is an addition doman before unrecognized HisKA: d1, <hole>, <HisKA>, HTPAse_c	
	if HisKAspace >= (minLengthForHisKA + minDomainLength):
		addHoleAndDomain(domainsOutput, domainPenultimate["env_to"]+1, holeEnd, HisKAdomain)
		domainsOutput.append(domainUltimate)
	else:
		domainsOutput.extend([HisKAdomain, domainUltimate])

def checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength, coord=1):
	firstDomain = domains[0]
	#process first domain
	if firstDomain["env_from"] > minDomainLength:
		addHoleAndDomain(domainsOutput, coord, firstDomain["env_from"]-1, firstDomain)
	else:
		domainsOutput.append(firstDomain)
	#process the rest of the domains
	for domain in domains[1:]:
		if (domain["env_from"] - firstDomain["env_to"]) >= minDomainLength:
			addHoleAndDomain(domainsOutput, firstDomain["env_to"]+1, domain["env_from"]-1, domain)
		else:
			domainsOutput.append(domain)
		firstDomain = domain

def addHoleAndDomain(domainsOutput, holeStart, holeEnd, domain):
	domainsOutput.extend([{"name": "hole", "env_from": holeStart, "env_to": holeEnd}, domain])

#### Process holes and domains END ####

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
	processDomains()

main(sys.argv)
