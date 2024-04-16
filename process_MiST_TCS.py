#!/usr/bin/python3
import sys, getopt
import urllib.request, urllib.parse, urllib.error
import json
import collections
import os.path

USAGE = "\n\nThe script queries MiST db via it's API for histidine kinases in genomes.\n\n" + \
	" \n\n" + \
	"python 	" + sys.argv[0] + '''
	-h || --help               - help
	-i || --ifile              - input file
	-o || --ofile              - output file
'''

GENOMES_URL = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes"
SIGNAL_GENES_HK = "/signal-genes?where.ranks=tcp,hk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"
SIGNAL_GENES_HHK = "/signal-genes?where.ranks=tcp,hhk&count&page=%PAGE%&per_page=100&fields.Gene.Aseq=pfam31"


INPUT_FILE = None
OUTPUT_FILE1 = None
OUTPUT_FILE2 = None

GENOME_VERSIONS = ""
TIMEOUT_FILE = "timeout_genomes.txt"

def initialize(argv):
	global INPUT_FILE, OUTPUT_FILE
	try:
		opts, args = getopt.getopt(argv[1:],"hi:s:o:p:t:",["help", "ifile=", "ofile="])
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

def initialyzeGenomeList():
	global GENOME_VERSIONS
	with open(INPUT_FILE, "r") as inputFile:
		GENOME_VERSIONS = [record.strip() for record in inputFile]

def retrieveSignalGenesFromMist(genomeVersion, additionaFieldsTemplate):
	genomeURL = GENOMES_URL + genomeVersion
	signalGeneList = list()
	noDataAnymore = False
	for num in range(1, 101):
		for iteration in range(0, 10):
			try:
				additionaFields = additionaFieldsTemplate.replace("%PAGE%", str(num))
				constructedUrl = genomeURL + additionaFields
				result = urllib.request.urlopen(constructedUrl)
				resultAsJson = json.loads(result.read().decode("utf-8"))
				if len(resultAsJson) == 0:   #No data anymore from this page on
					noDataAnymore = True
					break
				if "name" in resultAsJson:	#404 NotFoundError
					break
			#except ValueError: #504 Gateway timeouts
			except urllib.error.HTTPError:
			#except json.decoder.JSONDecodeError:   #504 Gateway timeouts  From Python 3.5+
				if iteration == 9:
					with open (TIMEOUT_FILE, "a") as timeoutFile:
						timeoutFile.write(genomeVersion + "\n")
				continue
			signalGeneList.extend(resultAsJson)
			break
		if noDataAnymore:
			break
	return signalGeneList			

def processDomains():
	genomeNumber = 1
	for genomeVersion in GENOME_VERSIONS:
		if genomeVersion.strip():
			print("Genome Number: " + str(genomeNumber))
			genomeNumber+=1
			listOfSignalGeneLists = []
			listOfSignalGeneLists.append(retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_HK))
			listOfSignalGeneLists.append(retrieveSignalGenesFromMist(genomeVersion, SIGNAL_GENES_HHK))
			domainCombinToCount = collections.defaultdict(int)

			for signalGeneList in listOfSignalGeneLists:
				for gene in signalGeneList:
					prepareDomains(gene, genomeVersion.strip(), domainCombinToCount)
			
			#domainCombinToCount is populated by this step and we can save the result
			domainCombinToCountList = sorted(list(domainCombinToCount.items()), key=lambda a: a[1], reverse=True)
			with open(OUTPUT_FILE2, "a") as outputFile:
				for domainCombinAndCount in domainCombinToCountList:
					outputFile.write(genomeVersion + "\t" + domainCombinAndCount[0] + "\t" + str(domainCombinAndCount[1]) + "\n")		

##*********************************************************************##
##********************** Domains processing block**********************##
def prepareDomains(gene, genomeVersion, domainCombinToCount):
	if "Gene" in gene and "Aseq" in gene["Gene"] and "pfam31" in gene["Gene"]["Aseq"]:
		#Ordering domains according to how they encoded in the gene
		domainsSorted = sorted(gene["Gene"]["Aseq"]["pfam31"], key=lambda x: x["ali_from"], reverse=False)
		if len (domainsSorted) > 0:
			domainsFiltered = removeOverlapps(domainsSorted)
			
			refseqVersion = gene["Gene"]["version"]
			geneStableId = gene["Gene"]["stable_id"]
			proteinLength = int(gene["Gene"]["length"]/3) - 1
			
			domainsOutput = processHoles(domainsFiltered)

			domainArchitecture = ""
			for domain in domainsOutput:
				domainArchitecture = domainArchitecture + "{}:{}-{}, ".format(domain["name"], domain["env_from"], domain["env_to"])
			
			#Generate a set of unique domain names and domain to count uniformly sorted
			#{'domain1': 1, 'domain2': 2}
			domainToCount = collections.defaultdict(int)
			for domain in domainsOutput:
				domainToCount[domain]+=1
			sortedDomainNames = sorted(domainToCount.keys())
			domainsFilteredNamesUniqueStr = ",".join(sortedDomainNames)
			domainsFilteredNamesUniqueCountsStr = ",".join(["{}:{}".format(domain, domainToCount[domain]) for domain in sortedDomainNames])

			#Save everything to file:
			with open(OUTPUT_FILE1, "a") as outputFile:
				outputRecord = "\t".join([genomeVersion, refseqVersion, geneStableId, proteinLength, domainArchitecture.rstrip(","), domainsFilteredNamesUniqueCountsStr, domainsFilteredNamesUniqueStr])
				outputFile.write(outputRecord + "\n")

			#Count a particular domain combination
			domainCombinToCount[domainsFilteredNamesUniqueStr]+=1
			
			return

	return False

def processHoles(domains):
	minDomainLength = 100
	minLengthForHisKA = 150
	domainsOutput = []
	if domains[-1] == 'HATPase_c':
		if domains[-2] == 'HisKA':
			checkAndAddHoles(domains, domainsOutput, minDomainLength)
		elif domains[-2] != 'HisKA':
			checkAndAddHoles(domains, domainsOutput, minDomainLength)
			#process the hole between the penultimate domain and the last domain
			HisKAspace = domains[-1]["env_from"] - domains[-2]["env_to"]
			#if the below condition is satisfied then there is an addition doman before HisKA	
			if HisKAspace >= (minLengthForHisKA + minDomainLength):
				holeEnd = domains[-2]["env_to"] + (HisKAspace - minLengthForHisKA) 
				addHoleAndDomain(domainsOutput, domains[-2]["env_to"]+1, holeEnd, domains[-1], True)
	elif domains[-1] != 'HATPase_c':
		#This means that the HATPase_c domain simply has not been recognized
		if domains[-1] == 'HisKA':
			checkAndAddHoles(domains, domainsOutput, minDomainLength)
		#This should not happen. Output the result if this happend:
		elif domains[-2] == 'HisKA':
			print ("Smething is not OK")
		#This should not happen either. Output the result if this happend:
		elif domains[-2] == 'HATPase_c':
			print ("Smething is not OK")
		#This should not happen at all:
		else:
			print ("Smething is completely wrong")
	return domainsOutput

def checkAndAddHoles(domains, domainsOutput, minDomainLength):
	firstDomain = domains[0]
	isHisKAorHATPase = False
	#process first domain
	if firstDomain["env_from"] > minDomainLength:
		addHoleAndDomain(domainsOutput, 1, firstDomain["env_from"]-1, firstDomain, isHisKAorHATPase)
	#process middle domains; two last domains are HisKA and HATPase_c and do not require special processing
	for domain in domains[1:-1]:
		if domain["name"] == "HisKA":
			isHisKAorHATPase = True
		if (domain["env_from"] - firstDomain["env_to"]) >= minDomainLength:
			addHoleAndDomain(domainsOutput, firstDomain["env_to"]+1, domain["env_from"]-1, domain, isHisKAorHATPase)
			firstDomain = domain	

def addHoleAndDomain(domainsOutput, holeStart, holeEnd, domain, isHisKAorHATPase=False):
	domainsOutput.append({"name": "hole", "env_from": holeStart, "env_to": holeEnd})
	if not isHisKAorHATPase:
		domainsOutput.append(domain)


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
	initialyzeGenomeList()

main(sys.argv)
	
