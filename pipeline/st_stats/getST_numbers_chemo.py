#!/usr/bin/python3
from collections import OrderedDict
import sys
import getopt

USAGE = "\n\nThe script summarizes statistic of chemotacits systems created by the process_MiST_ST_counts.py script.\n\n" + \
    "python 	" + sys.argv[0] + '''
    -h || --help               - help
    -a || --afile              - input file (CheA statistics file). Content example:
            genomeID	$total	chew	checx	other	F1	F2	F3	F4	F5	F6	F7	F8	F9	F10	F11	F12	F13	F14	F15	F16	F17	Acf	Tfp	Taxonomy
            GCA_001872605.1	1	1	2	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	d__Bacteria;p__Armatimonadota;c__Abditibacteria;o__CG2-30-59-28;f__CG2-30-59-28;g__CG2-30-59-28;s__CG2-30-59-28 sp001872605
            GCA_001873295.1	2	2	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	d__Bacteria;p__CG2-30-70-394;c__CG2-30-70-394;o__CG2-30-70-394;f__CG2-30-70-394;g__CG2-30-70-394;s__CG2-30-70-394 sp001873295
            GCA_002717565.1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	d__Bacteria;p__Chloroflexota;c__Dehalococcoidia;o__GCA-2717565;f__GCA-2717565;g__GCA-2717565;s__GCA-2717565 sp002717565
    -b || --bfile              - input file (CheB statistics file)
    -r || --rfile              - input file (CheR statistics file) 
    -z || --zfile              - input file (CheZ statistics file)
    -d || --dfile              - input file (CheD statistics file)
    -v || --vfile              - input file (CheV statistics file)
    -m || --mfile              - input file (MCP statistics file)
    -t || --taxlevel           - taxonomy level for summarization. One of: species, genus, family, order, class, taxlevel, kingdom, or acorss. 'across' meas across all phyla
    -o || --ofile              - output file
    '''

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
INPUT_FILE_TO_COMPONENT = {CHEA_FILE_INPUT: "CheA", CHEB_FILE_INPUT: "CheB", CHER_FILE_INPUT: "CheR", CHEZ_FILE_INPUT: "CheZ", CHED_FILE_INPUT: "CheD", CHEV_FILE_INPUT: "CheV", MCP_FILE_INPUT: "MCP"}


GENOME_TO_CHEA_NUMBER = {}

TAXLEVEL = "phylum"
TAXONOMY_TO_LEVEL = {"species": 7, "genus": 6, "family": 5, "order": 4, "class": 3, "phylum": 2, "kingdom": 1}

#TAXONOMY_LEVEL_TO_COMPONENT_TO_DATA = {"pylum1": {"CheA": {"total": 0, "F1": 12, ...}, "CheB": {"total": 0, "F1": 12, ...}, }
TAXONOMY_LEVEL_TO_COMPONENT_TO_DATA = {}

def initialize(argv):
    global A_FILE, B_FILE, R_FILE, Z_FILE, D_FILE, V_FILE, M_FILE, OUTPUT_FILE, TAXLEVEL
    try:
        opts, args = getopt.getopt(argv[1:],"hi:o:t:a:b:r:z:d:v:m:",["help", "afile=", "bfile=", "rfile=", "zfile=", "dfile=", "vfile=", "mfile=", "ofile=", "taxlevel="])
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
            elif opt in ("-a", "--afile"):
                A_FILE = str(arg).strip()
            elif opt in ("-b", "--bfile"):
                B_FILE = str(arg).strip()
            elif opt in ("-r", "--rfile"):
                R_FILE = str(arg).strip()
            elif opt in ("-z", "--zfile"):
                Z_FILE = str(arg).strip()
            elif opt in ("-d", "--dfile"):
                D_FILE = str(arg).strip()
            elif opt in ("-v", "--vfile"):
                V_FILE = str(arg).strip()
            elif opt in ("-m", "--mfile"):
                M_FILE = str(arg).strip()
            elif opt in ("-o", "--ofile"):
                OUTPUT_FILE = str(arg).strip()
            elif opt in ("-t", "--taxlevel"):
                TAXLEVEL = str(arg).strip().lower()
                OUTPUT_FILE += "_" + TAXLEVEL
    except Exception as e:
        print(("===========ERROR==========\n " + str(e) + USAGE))
        sys.exit(2)

def processSTstatistics(fileToProcess, data):
    #taxlevel_to_data = {"taxlevel1": {"total": 0, ...}}
    taxlevel_to_data = {}
    taxlevel = "Across_Phyla"
    with open (fileToProcess, "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            line = line.strip().split("\t")
            if lineNumber > 0:
                if fileToProcess == CHEA_FILE_INPUT:
                    GENOME_TO_CHEA_NUMBER[line[0]] = int(line[1])
                if GENOME_TO_CHEA_NUMBER[line[0]] >= 1:
                    if TAXLEVEL in TAXONOMY_TO_LEVEL:
                        taxlevel = ";".join(line[-1].split(";")[:TAXONOMY_TO_LEVEL[TAXLEVEL]])
                    else:
                        taxlevel = "Across_Phyla"
                    if lineNumber == 1 or (taxlevel not in taxlevel_to_data and lineNumber > 1):
                        taxlevel_to_data[taxlevel] = data.copy()

                #I use num+1 becuase the first column in the input file is a genome identifier,
                #while elements of INPUT_FILE_TO_DATA[fileToProcess] look like: "total", "F1, F2", etc
                #There we add 1 to shift by one position after genome identifier. Look at the input file format in USAGE
                for num, field in enumerate(INPUT_FILE_TO_DATA[fileToProcess]):
                    TAXONOMY_LEVEL_TO_COMPONENT_TO_DATA[taxlevel][INPUT_FILE_TO_COMPONENT[fileToProcess]]["recordNumber"] += 1
                    TAXONOMY_LEVEL_TO_COMPONENT_TO_DATA[taxlevel][INPUT_FILE_TO_COMPONENT[fileToProcess]][field] += int(line[num+1])

    return taxlevel_to_data
                    
def finalizeDataAndPrint(taxlevel_to_data, inputFile):
    with open (inputFile+"_summary.txt", "w") as outputFile:
        for taxlevel, data in taxlevel_to_data.items():
            for param in data:
                if param != "recordNumber":
                    data[param] = data[param]/data["recordNumber"]
            dataValues = map(roundToFirstDecim, data.values())
            outputFile.write(taxlevel + "\t" + "\t".join(map(str, dataValues)) + "\n")

def roundToFirstDecim(value):
    return round(value, 1)
	
def main(argv):
    initialize(argv)
    for inputFile, data in INPUT_FILE_TO_DATA.items():
        print ("axLevel" + "\t" + "\t".join(data.keys()))
        taxlevel_to_data = processSTstatistics(inputFile, data)
        finalizeDataAndPrint(taxlevel_to_data, inputFile)

main(sys.argv)