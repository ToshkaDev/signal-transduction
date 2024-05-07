#!/usr/bin/python3
import sys, getopt
import urllib.request, urllib.parse, urllib.error
import json
import collections
import os.path
import time

USAGE = "\n\nThe script queries MiST db via it's API for histidine kinases in genomes.\n\n" + \
	" \n\n" + \
	"python 	" + sys.argv[0] + '''
	-h || --help               - help
	-i || --ifile              - input file
	-f || --ffile              - first output file
	-s || --sfile              - second output file
'''

GENOMES_URL = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/"
STP_MATRIX = "/stp-matrix"
SIGNAL_GENES_HK = "/signal-genes?where.component_id=%COMPONENT_ID%&where.ranks=tcp,hk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"
SIGNAL_GENES_HHK = "/signal-genes?where.component_id=%COMPONENT_ID%&where.ranks=tcp,hhk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"
#SIGNAL_GENES_HK = "/signal-genes?where.ranks=tcp,hk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"
#SIGNAL_GENES_HHK = "/signal-genes?where.ranks=tcp,hhk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"


INPUT_FILE = None
OUTPUT_FILE1 = None
OUTPUT_FILE2 = None

GENOME_VERSIONS = ""
TIMEOUT_FILE = "timeout_genomes.txt"

HIS_KINASE_DIM_DOMAINS = ["HisKA", "HisKA_2", "HisKA_3", "H-kinase_dim", "His_kinase"]
HIS_KINASE_CATAL_DOMAINS = ["HATPase_c", "HATPase_c_2", "HATPase_c_5", "HWE_HK"]

def initialize(argv):
	global INPUT_FILE, OUTPUT_FILE1, OUTPUT_FILE2, GENOME_VERSIONS
	try:
		opts, args = getopt.getopt(argv[1:],"hi:f:s:",["help", "ifile=", "ffile=", "sfile="])
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
	except Exception as e:
		print("===========ERROR==========\n " + str(e) + USAGE)
		sys.exit(2)
	
	with open(INPUT_FILE, "r") as inputFile:
		GENOME_VERSIONS = [record.strip() for record in inputFile]

#The function first retreives stp-matrix and after that signal genes for those components that have target signaling systems
#This is done, because it is much faster this way
def retrieveSignalGenesFromMist(genomeVersion, additionaFieldsTemplate):
	if additionaFieldsTemplate == SIGNAL_GENES_HK:
		hkOrhhK = "hk"
	elif additionaFieldsTemplate == SIGNAL_GENES_HHK:
		hkOrhhK = "hhk"
	genomeURL = GENOMES_URL + genomeVersion
	noDataAnymore = False

	#Get stp-matrix and look at the numnber of components and save those components that have two-component systems.
	#They will be saved in componentsWithTcp list.
	componentsWithTcp = list() #componentsWithTcp will be poppulated
	signalGenesRetriever(genomeURL + STP_MATRIX, componentsWithTcp, genomeVersion, hkOrhhK, noDataAnymore)
	#Retrieve signal genes in those components that have two-component systems
	signalGeneList = list()
	for component in componentsWithTcp:
		for num in range(1, 101):
			additionaFields = additionaFieldsTemplate.replace("%COMPONENT_ID%", str(component["id"])).replace("%PAGE%", str(num))
			constructedUrl = genomeURL + additionaFields
			noDataAnymore = signalGenesRetriever(constructedUrl, signalGeneList, genomeVersion, False, noDataAnymore)
			if noDataAnymore:
				break
	return signalGeneList			

def signalGenesRetriever(url, elementList, genomeVersion, tcpMatrix, noDataAnymore):
	for iteration in range (0, 10):
		try:
			result = urllib.request.urlopen(url)
			resultAsJson = json.loads(result.read().decode("utf-8"))
			if len(resultAsJson) == 0:   #No data anymore from this page on
				noDataAnymore = True
				break
			if "name" in resultAsJson:	#404 NotFoundError
				break
		except urllib.error.HTTPError:
		#except json.decoder.JSONDecodeError:   #504 Gateway timeouts  From Python 3.5+
			if iteration == 9:
				with open (TIMEOUT_FILE, "a") as timeoutFile:
					timeoutFile.write(genomeVersion + "\n")
			#sleep 5 seconds if gateway timeout happened		
			time.sleep(5)
			continue

		if tcpMatrix and "tcp" in resultAsJson["counts"]:
			for component in resultAsJson["components"]:
				if "tcp" in component["counts"] and tcpMatrix in component["counts"]["tcp"]:
					elementList.append(component)
		else:
			elementList.extend(resultAsJson)
		break

	return noDataAnymore

def processDomains():
	genomeNumber = 1
	for genomeVersion in GENOME_VERSIONS:
		genomeVersion = genomeVersion.strip()
		if genomeVersion:
			print("Genome Number: " + str(genomeNumber))
			genomeNumber+=1
			listOfSignalGeneLists = []
			listOfSignalGeneLists.append(retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_HK))
			listOfSignalGeneLists.append(retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_HHK))
			domainCombinToCount = collections.defaultdict(int)

			for signalGeneList in listOfSignalGeneLists:
				for gene in signalGeneList:
					prepareDomains(gene, genomeVersion, domainCombinToCount)
			
			#domainCombinToCount in each genome is populated at the previous step and we can save the result
			domainCombinToCountList = sorted(list(domainCombinToCount.items()), key=lambda a: a[1], reverse=True)
			with open(OUTPUT_FILE2, "a") as outputFile:
				for domainCombinAndCount in domainCombinToCountList:
					outputFile.write(genomeVersion + "\t" + domainCombinAndCount[0] + "\t" + str(domainCombinAndCount[1]) + "\n")		

##*********************************************************************##
##********************** Domains processing block**********************##
def prepareDomains(gene, genomeVersion, domainCombinToCount):
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
			print (geneStableId)
			domainsOutput = processHoles(domainsFiltered)

			domainArchitecture = ""
			domainArchitectureSensDomsOnly = ""
			for domain in domainsOutput:
				domainArchitecture = domainArchitecture + "{}:{}-{},".format(domain["name"], domain["env_from"], domain["env_to"])
				if domain["name"] != "hole" and domain["name"].lstrip("<").rstrip(">") not in HIS_KINASE_CATAL_DOMAINS and domain["name"].lstrip("<").rstrip(">") not in HIS_KINASE_DIM_DOMAINS:
					domainArchitectureSensDomsOnly = domainArchitectureSensDomsOnly + "," + domain["name"]

			#Generate a set of unique domain names and domain to count uniformly sorted
			#{'domain1': 1, 'domain2': 2}
			domainToCount = collections.defaultdict(int)
			for domain in domainsOutput:
				domainToCount[domain["name"]]+=1
			sortedDomainNames = sorted(domainToCount.keys())
			domainsFilteredNamesUniqueStr = ",".join(sortedDomainNames)
			domainsFilteredNamesUniqueCountsStr = ",".join(["{}:{}".format(domain, domainToCount[domain]) for domain in sortedDomainNames])

			#Save everything to file:
			with open(OUTPUT_FILE1, "a") as outputFile:
				outputRecord = "\t".join([genomeVersion, refseqVersion, geneStableId, proteinLength, domainArchitecture.rstrip(","), domainArchitectureSensDomsOnly.lstrip(","), domainsFilteredNamesUniqueCountsStr, domainsFilteredNamesUniqueStr])
				outputFile.write(outputRecord + "\n")

			#Count a particular domain combination
			domainCombinToCount[domainsFilteredNamesUniqueStr]+=1
			
			return

	return False

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
			#if not no domains upstream of the HATPase domain
			checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength)
	#if HisKA is present but HATPase is not (as HATPase cases are process above)
	elif HisKA_ind >=0:
		checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength)
	#if both HisKA and HATPase are not in my two lists, HIS_KINASE_DIM_DOMAINS and HIS_KINASE_CATAL_DOMAINS, still process all the domains:
	else:
		checkAndAddHolesAndDomains(domains, domainsOutput, minDomainLength)
		
	return domainsOutput

def HisKAprocessing(domainUltimate, domainPenultimate, minLengthForHisKA, minDomainLength, domainsOutput):
	#process the hole between the penultimate domain and the HATPase domain, under the assumption that between these two domains the HisKA domain can be present
	#      domainPenultimate       domainUltimate 
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
	
