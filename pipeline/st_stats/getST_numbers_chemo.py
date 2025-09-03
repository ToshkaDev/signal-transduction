#!/usr/bin/python3
from collections import OrderedDict

TEMPLATE_TO_COUNTS = OrderedDict([("$total", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), ("recordNumber", 0)])
CHEA_TO_COUNTS = OrderedDict([("$total", 0), ("chew", 0), ("checx", 0), ("other", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), ("recordNumber", 0)])
CHED_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEZ_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEV_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEB_TO_COUNTS = OrderedDict([("$total", 0), ("MAC1", 0), ("MAC2", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), ("recordNumber", 0)])
CHER_TO_COUNTS = CHEB_TO_COUNTS.copy()
MCP_TO_COUNTS = OrderedDict([("$total", 0), ("64H", 0), ("58H", 0), ("52H", 0), ("48H", 0), ("44H", 0), ("42H", 0), ("40H", 0), ("38H", 0), ("36H", 0), ("34H", 0), ("28H", 0), ("24H", 0), ("recordNumber", 0)]) 

CHEA_FILE_INPUT = "CheA_across_genomes_fullTaxonomy_shortToPhylum.txt"
CHEB_FILE_INPUT = "CheB_across_genomes_fullTaxonomy_shortToPhylum.txt"
CHER_FILE_INPUT = "CheR_across_genomes_fullTaxonomy_shortToPhylum.txt"
CHEZ_FILE_INPUT = "CheZ_across_genomes_fullTaxonomy_shortToPhylum.txt"
CHED_FILE_INPUT = "CheD_across_genomes_fullTaxonomy_shortToPhylum.txt"
CHEV_FILE_INPUT = "CheV_across_genomes_fullTaxonomy_shortToPhylum.txt"
MCP_FILE_INPUT = "MCP_across_genomes_fullTaxonomy_shortToPhylum.txt"

INPUT_FILE_TO_DATA = OrderedDict([(CHEA_FILE_INPUT, CHEA_TO_COUNTS), (CHEB_FILE_INPUT, CHEB_TO_COUNTS), (CHER_FILE_INPUT, CHER_TO_COUNTS), (CHEZ_FILE_INPUT, CHEZ_TO_COUNTS), (CHED_FILE_INPUT, CHED_TO_COUNTS), (CHEV_FILE_INPUT, CHEV_TO_COUNTS), (MCP_FILE_INPUT, MCP_TO_COUNTS)])

GENOME_TO_CHEA_NUMBER = {}
#GenomeId(0)    Total(1)   OCP_total(2)    TCP_total(3) TCP_HK(4)   TCP_HHK(5)  TCP_RR(6)   TCP_HRR(7)  TCP_Other(8)    ChemSys(9) ECF(10) Other(11)   Taxonomy-Phylum(12)
def processSTstatistics(fileToProcess, data):
    #phylum_to_data = {"pylum1": {"total": 0, ...}}
    phylum_to_data = {}
    phylum = "Across_Phyla"
    with open (fileToProcess, "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            line = line.strip().split("\t")
            if lineNumber > 0:
                if fileToProcess == CHEA_FILE_INPUT:
                    GENOME_TO_CHEA_NUMBER[line[0]] = int(line[1])
                if GENOME_TO_CHEA_NUMBER[line[0]] >= 1:
                    if "d__" in line[-1]:
                        phylum = line[-1]
                    else:
                        phylum = "Across_Phyla"   
                    if lineNumber == 1 or (phylum not in phylum_to_data and lineNumber > 1):
                        phylum_to_data[phylum] = data.copy()
                        

                    #increment statistics 
                    #I use num+1 becuase the first column in the input file is a genome identifier
                    for num, component in enumerate(phylum_to_data[phylum]):                     
                        if component != "recordNumber":
                            if fileToProcess != CHEA_FILE_INPUT:
                                #Devide the numbero of components by the number of chemotaxis systems based on the numbr of cheAs for a given genome
                                phylum_to_data[phylum][component] += float(line[num+1])/GENOME_TO_CHEA_NUMBER[line[0]]
                            else:
                                phylum_to_data[phylum][component] += float(line[num+1])
                        else:
                            phylum_to_data[phylum][component]+=1
          
    return phylum_to_data
                    
def finalizeDataAndPrint(phylum_to_data, inputFile):
    with open (inputFile+"_summary.txt", "w") as outputFile:
        for phylum, data in phylum_to_data.items():
            for param in data:
                if param != "recordNumber":
                    data[param] = data[param]/data["recordNumber"]
            dataValues = map(roundToFirstDecim, data.values())
            outputFile.write(phylum + "\t" + "\t".join(map(str, dataValues)) + "\n")

def roundToFirstDecim(value):
    return round(value, 1)
	
def main():
    for inputFile, data in INPUT_FILE_TO_DATA.items():
        print ("Phylum" + "\t" + "\t".join(data.keys()))
        phylum_to_data = processSTstatistics(inputFile, data)
        finalizeDataAndPrint(phylum_to_data, inputFile)

main()