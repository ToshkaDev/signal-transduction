#!/usr/bin/python3
from collections import OrderedDict
import sys
import getopt

USAGE = "\n\nThe script summarizes statistic of chemotaxis systems created by the process_MiST_ST_counts.py script.\n" + \
    "For each taxonomy level (species, genus, family, order, class, phylum, kingdom) it calculates the average number of each chemotaxis system component per genome.\n" + \
    "python 	" + sys.argv[0] + '''
    -h || --help               - help
    -a || --afile              - input file (CheA statistics file).    
    -b || --bfile              - input file (CheB statistics file)
    -r || --rfile              - input file (CheR statistics file) 
    -z || --zfile              - input file (CheZ statistics file)
    -d || --dfile              - input file (CheD statistics file)
    -v || --vfile              - input file (CheV statistics file)
    -m || --mfile              - input file (MCP statistics file)
    -t || --taxlevel           - taxonomy level for summarization. One of: species, genus, family, order, class, taxlevel, kingdom, or acorss. 'across' means across all phyla.

    Input files are tab delimited and have the following format (ex., -a file, CheA statistics file):
    genomeID	$total	chew	checx	other	F1	F2	F3	F4	F5	F6	F7	F8	F9	F10	F11	F12	F13	F14	F15	F16	F17	Acf	Tfp	Taxonomy
    GCA_001872605.1	1	1	2	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	d__Bacteria;p__Armatimonadota;c__Abditibacteria;o__CG2-30-59-28;f__CG2-30-59-28;g__CG2-30-59-28;s__CG2-30-59-28 sp001872605
    GCA_001873295.1	2	2	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	d__Bacteria;p__CG2-30-70-394;c__CG2-30-70-394;o__CG2-30-70-394;f__CG2-30-70-394;g__CG2-30-70-394;s__CG2-30-70-394 sp001873295
    GCA_002717565.1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	d__Bacteria;p__Chloroflexota;c__Dehalococcoidia;o__GCA-2717565;f__GCA-2717565;g__GCA-2717565;s__GCA-2717565 sp002717565
    '''

TEMPLATE_TO_COUNTS = OrderedDict([("$total", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), ("recordNumber", 0)])
CHEA_TO_COUNTS = OrderedDict([("$total", 0), ("chew", 0), ("checx", 0), ("other", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), ("recordNumber", 0)])
CHED_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEZ_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEV_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEB_TO_COUNTS = OrderedDict([("$total", 0), ("MAC1", 0), ("MAC2", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), ("recordNumber", 0)])
CHER_TO_COUNTS = CHEB_TO_COUNTS.copy()
MCP_TO_COUNTS = OrderedDict([("$total", 0), ("64H", 0), ("58H", 0), ("52H", 0), ("48H", 0), ("44H", 0), ("42H", 0), ("40H", 0), ("38H", 0), ("36H", 0), ("34H", 0), ("28H", 0), ("24H", 0), ("recordNumber", 0)]) 

A_FILE = "process_MiST_ST_counts_CheA_across_genomes.txt"
B_FILE = "process_MiST_ST_counts_CheB_across_genomes.txt"
R_FILE = "process_MiST_ST_counts_CheR_across_genomes.txt"
Z_FILE = "process_MiST_ST_counts_CheZ_across_genomes.txt"
D_FILE = "process_MiST_ST_counts_CheD_across_genomes.txt"
V_FILE = "process_MiST_ST_counts_CheV_across_genomes.txt"
MCP_FILE = "process_MiST_ST_counts_MCP_across_genomes.txt"

INPUT_FILE_TO_DATA = OrderedDict([(A_FILE, CHEA_TO_COUNTS), (B_FILE, CHEB_TO_COUNTS), (R_FILE, CHER_TO_COUNTS), (Z_FILE, CHEZ_TO_COUNTS), (D_FILE, CHED_TO_COUNTS), (V_FILE, CHEV_TO_COUNTS), (MCP_FILE, MCP_TO_COUNTS)])

GENOME_TO_CHEA_NUMBER = {}

TAXLEVEL = "phylum"
TAXONOMY_TO_LEVEL = {"species": 7, "genus": 6, "family": 5, "order": 4, "class": 3, "phylum": 2, "kingdom": 1}

def initialize(argv):
    global A_FILE, B_FILE, R_FILE, Z_FILE, D_FILE, V_FILE, MCP_FILE, TAXLEVEL
    try:
        opts, args = getopt.getopt(argv[1:],"hi:t:a:b:r:z:d:v:m:",["help", "afile=", "bfile=", "rfile=", "zfile=", "dfile=", "vfile=", "mfile=", "taxlevel="])
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
                MCP_FILE = str(arg).strip()
            elif opt in ("-t", "--taxlevel"):
                TAXLEVEL = str(arg).strip().lower()
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
                if TAXLEVEL in TAXONOMY_TO_LEVEL:
                    taxlevel = ";".join(line[-1].split(";")[:TAXONOMY_TO_LEVEL[TAXLEVEL]])
                else:
                    taxlevel = "Across_Phyla"
                if lineNumber == 1 or (taxlevel not in taxlevel_to_data and lineNumber > 1):
                    taxlevel_to_data[taxlevel] = data.copy()

                #I use num+1 becuase the first column in the input file is a genome identifier,
                #while elements of INPUT_FILE_TO_DATA[fileToProcess] look like: "total", "F1, F2", etc
                #There we add 1 to shift by one position after genome identifier. Look at the input file format in USAGE
                for num, field in enumerate(taxlevel_to_data[taxlevel]):
                    #-2 because last column is taxonomy, first is genomeID
                    if num < len(line)-2:  
                        taxlevel_to_data[taxlevel]["recordNumber"] += 1
                        taxlevel_to_data[taxlevel][field] += int(line[num+1])

    return taxlevel_to_data

def finalizeDataAndSave(taxlevel_to_data, inputFile, data):
    with open (inputFile.split(".")[0]+"_"+TAXLEVEL, "a") as output_file:
        output_file.write("taxLevel" + "\t" + "\t".join(data.keys()) + "\n")

        for taxlevel, data in taxlevel_to_data.items():
            for param in data:
                if param != "recordNumber":
                    data[param] = data[param]/data["recordNumber"]
            dataValues = map(roundToFirstDecim, data.values())
            output_file.write(taxlevel + "\t" + "\t".join(map(str, dataValues)) + "\n")

def roundToFirstDecim(value):
    return round(value, 1)
	
def main(argv):
    initialize(argv)
    for inputFile, data in INPUT_FILE_TO_DATA.items():
        print ("taxLevel" + "\t" + "\t".join(data.keys()))
        taxlevel_to_data = processSTstatistics(inputFile, data)
        finalizeDataAndSave(taxlevel_to_data, inputFile, data)

main(sys.argv)