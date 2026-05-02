#!/usr/bin/python3
from collections import OrderedDict
import sys
import getopt
import re

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
    -s || --sfile              - input file (GTDB taxonomy metadata file)
    -t || --taxlevel           - taxonomy level for summarization. One of: species, genus, family, order, class, taxlevel, kingdom, or acorss. 'across' means across all phyla.

    Input files are tab delimited and have the following format (ex., -a file, CheA statistics file):
    genome	genome_accession	total	chew	checx	other	F1	F2	F3	F4	F5	F6	F7	F8	F9	F10	F11	F12	F13	F14	F15	F16	F17	Acf	Tfp
    GCA_001872605.1	GCA_001872605	1	1	2	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    GCA_001873295.1	GCA_001873295	2	2	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    '''

TEMPLATE_TO_COUNTS = OrderedDict([("total", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), \
                                  	("total_norm_by_protein_count", 0), ("F1_norm_by_protein_count", 0), ("F2_norm_by_protein_count", 0), ("F3_norm_by_protein_count", 0), ("F4_norm_by_protein_count", 0), ("F5_norm_by_protein_count", 0), \
									("F6_norm_by_protein_count", 0), ("F7_norm_by_protein_count", 0), ("F8_norm_by_protein_count", 0), ("F9_norm_by_protein_count", 0), ("F10_norm_by_protein_count", 0), ("F11_norm_by_protein_count", 0), \
									("F12_norm_by_protein_count", 0), ("F13_norm_by_protein_count", 0), ("F14_norm_by_protein_count", 0), ("F15_norm_by_protein_count", 0), ("F16_norm_by_protein_count", 0), ("F17_norm_by_protein_count", 0), \
									("Acf_norm_by_protein_count", 0), ("Tfp_norm_by_protein_count", 0), ("recordNumber", 0)])
CHEA_TO_COUNTS = OrderedDict([("total", 0), ("chew", 0), ("checx", 0), ("other", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), \
                              	("total_norm_by_protein_count", 0), ("chew_norm_by_protein_count", 0), ("checx_norm_by_protein_count", 0), ("other_norm_by_protein_count", 0), ("F1_norm_by_protein_count", 0), ("F2_norm_by_protein_count", 0), \
								("F3_norm_by_protein_count", 0), ("F4_norm_by_protein_count", 0), ("F5_norm_by_protein_count", 0), ("F6_norm_by_protein_count", 0), ("F7_norm_by_protein_count", 0), ("F8_norm_by_protein_count", 0), ("F9_norm_by_protein_count", 0), \
								("F10_norm_by_protein_count", 0), ("F11_norm_by_protein_count", 0), ("F12_norm_by_protein_count", 0), ("F13_norm_by_protein_count", 0), ("F14_norm_by_protein_count", 0), ("F15_norm_by_protein_count", 0), \
								("F16_norm_by_protein_count", 0), ("F17_norm_by_protein_count", 0), ("Acf_norm_by_protein_count", 0), ("Tfp_norm_by_protein_count", 0), ("recordNumber", 0)])
CHED_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEZ_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEV_TO_COUNTS = TEMPLATE_TO_COUNTS.copy()
CHEB_TO_COUNTS = OrderedDict([("total", 0), ("MAC1", 0), ("MAC2", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("Acf", 0), ("Tfp", 0), \
                              	("total_norm_by_protein_count", 0), ("MAC1_norm_by_protein_count", 0), ("MAC2_norm_by_protein_count", 0), ("F1_norm_by_protein_count", 0), ("F2_norm_by_protein_count", 0), ("F3_norm_by_protein_count", 0), \
								("F4_norm_by_protein_count", 0), ("F5_norm_by_protein_count", 0), ("F6_norm_by_protein_count", 0), ("F7_norm_by_protein_count", 0), ("F8_norm_by_protein_count", 0), ("F9_norm_by_protein_count", 0), ("F10_norm_by_protein_count", 0), \
								("F11_norm_by_protein_count", 0), ("F12_norm_by_protein_count", 0), ("F13_norm_by_protein_count", 0), ("F14_norm_by_protein_count", 0), ("F15_norm_by_protein_count", 0), ("F16_norm_by_protein_count", 0), \
								("F17_norm_by_protein_count", 0), ("Acf_norm_by_protein_count", 0), ("Tfp_norm_by_protein_count", 0), ("recordNumber", 0)])
CHER_TO_COUNTS = CHEB_TO_COUNTS.copy()
MCP_TO_COUNTS = OrderedDict([("total", 0), ("64H", 0), ("58H", 0), ("52H", 0), ("48H", 0), ("44H", 0), ("42H", 0), ("40H", 0), ("38H", 0), ("36H", 0), ("34H", 0), ("28H", 0), ("24H", 0), \
                             	("total_norm_by_protein_count", 0), ("64H_norm_by_protein_count", 0), ("58H_norm_by_protein_count", 0), ("52H_norm_by_protein_count", 0), \
								("48H_norm_by_protein_count", 0), ("44H_norm_by_protein_count", 0), ("42H_norm_by_protein_count", 0), ("40H_norm_by_protein_count", 0), \
								("38H_norm_by_protein_count", 0), ("36H_norm_by_protein_count", 0), ("34H_norm_by_protein_count", 0), ("28H_norm_by_protein_count", 0), ("24H_norm_by_protein_count", 0), ("recordNumber", 0)]) 

GTDB_FILE = None
A_FILE = "get_counts_CheA_across_genomes.txt"
B_FILE = "get_counts_CheB_across_genomes.txt"
R_FILE = "get_counts_CheR_across_genomes.txt"
Z_FILE = "get_counts_CheZ_across_genomes.txt"
D_FILE = "get_counts_CheD_across_genomes.txt"
V_FILE = "get_counts_CheV_across_genomes.txt"
MCP_FILE = "get_counts_MCP_across_genomes.txt"

INPUT_FILE_TO_DATA = OrderedDict([(A_FILE, CHEA_TO_COUNTS), (B_FILE, CHEB_TO_COUNTS), (R_FILE, CHER_TO_COUNTS), (Z_FILE, CHEZ_TO_COUNTS), (D_FILE, CHED_TO_COUNTS), (V_FILE, CHEV_TO_COUNTS), (MCP_FILE, MCP_TO_COUNTS)])

GENOME_TO_CHEA_NUMBER = {}

TAXLEVEL = "phylum"
TAXONOMY_TO_LEVEL = {"species": 7, "genus": 6, "family": 5, "order": 4, "class": 3, "phylum": 2, "kingdom": 1}
#Major mode ("mm") table (see the MiST database for details)
PROTEIN_TYPE = "CheA"
INPUT_FILE_TO_PROTEIN_TYPE = dict([(A_FILE, "CheA"), (B_FILE, "CheB"), (R_FILE, "CheR"), (Z_FILE, "CheZ"), (D_FILE, "CheD"), (V_FILE, "CheV"), (MCP_FILE, "MCP")])
GENOME_TO_TAXONOMY = {}
ADDITIONAL_HEADERS = "\t".join(["gtdb_taxonomy_string", "gtdb_taxonomy_last", "gtdb_taxonomy_rank", "protein_type"])

def initialize(argv):
    global INPUT_FILE_TO_DATA, GTDB_FILE, A_FILE, B_FILE, R_FILE, Z_FILE, D_FILE, V_FILE, MCP_FILE, TAXLEVEL, PROTEIN_TYPE
    try:
        opts, args = getopt.getopt(argv[1:],"hi:s:a:b:r:z:d:v:m:t:",["help", "sfile=", "afile=", "bfile=", "rfile=", "zfile=", "dfile=", "vfile=", "mfile=", "taxlevel="])
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
            elif opt in ("-s", "--sfile"):
                GTDB_FILE = str(arg).strip()
            elif opt in ("-t", "--taxlevel"):
                TAXLEVEL = str(arg).strip().lower()
    except Exception as e:
        print(("===========ERROR==========\n " + str(e) + USAGE))
        sys.exit(2)
    
    # regex to remove prefix .__ (like d__, p__, etc. from the GTDB taxonomomy string)
    regex=r'.__'
    #taxlevel_to_data = {"taxlevel1": {"total": 0, ...}}
    with open(GTDB_FILE, "r") as iFile2:
        for line in iFile2:
            #$0 - genome version, $1 - genome accession, $2 - genome size, $3 - protein counts, $4 - GTDB taxonomy, $5 - NCBI taxonomy
            record = line.strip().split("\t")
            genome_version = record[0]
            taxonomy = ";".join(record[4].split(";")[:TAXONOMY_TO_LEVEL[TAXLEVEL]])
            taxonomy = re.sub(regex, "", taxonomy)
            GENOME_TO_TAXONOMY[genome_version] = taxonomy

    INPUT_FILE_TO_DATA = OrderedDict([(A_FILE, CHEA_TO_COUNTS), (B_FILE, CHEB_TO_COUNTS), (R_FILE, CHER_TO_COUNTS), (Z_FILE, CHEZ_TO_COUNTS), (D_FILE, CHED_TO_COUNTS), (V_FILE, CHEV_TO_COUNTS), (MCP_FILE, MCP_TO_COUNTS)])

def processSTstatistics(fileToProcess, data):
    taxlevel_to_data = {}
    with open (fileToProcess, "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            if lineNumber > 0:
                line = line.strip().split("\t")
                taxlevel = GENOME_TO_TAXONOMY[line[0]]

                if lineNumber == 1 or (lineNumber > 1 and taxlevel not in taxlevel_to_data):
                    taxlevel_to_data[taxlevel] = data.copy()

                taxlevel_to_data[taxlevel]["recordNumber"] += 1
                #I use num+2 becuase the first and 2d columns in the input file are a genome identifier and genome_accession,
                #while elements of INPUT_FILE_TO_DATA[fileToProcess] look like: "total", "F1, F2", etc
                #Therefore we add 2 to shift by two positions from the beginning of the line and land on the 3rd column. Look at the input file format in USAGE
                for idx, system in enumerate(taxlevel_to_data[taxlevel]):
                    taxlevel_to_data[taxlevel][system] += float(line[idx+2])

    return taxlevel_to_data

def finalizeDataAndSave(taxlevel_to_data, inputFile, data, protein_type):
    with open (inputFile.split(".")[0]+"_"+TAXLEVEL, "a") as output_file:
         # converting to list and write all the headers except for the last one, which is 'record_number'
        output_file.write(ADDITIONAL_HEADERS + "\t" + "\t".join(list(data.keys())[:-1]) + "\n")

        for taxlevel, data in taxlevel_to_data.items():
            for param in data:
                data[param] = data[param]/data["recordNumber"]
            dataValues = map(roundToFirstDecim, data.values())
            # Here converting the resulting map obejct to lost and writing all the values except for the last one (the number of records)
            output_file.write("\t".join([taxlevel, taxlevel.split(";")[-1], TAXLEVEL, protein_type])  + "\t" + "\t".join(list(map(str, dataValues))[:-1]) + "\n")

def roundToFirstDecim(value):
    return round(value, 1)
	
def main(argv):
    initialize(argv)
    for inputFile, data in INPUT_FILE_TO_DATA.items():
        print ("taxLevel" + "\t" + "\t".join(data.keys()))
        taxlevel_to_data = processSTstatistics(inputFile, data)
        finalizeDataAndSave(taxlevel_to_data, inputFile, data, INPUT_FILE_TO_PROTEIN_TYPE[inputFile])

main(sys.argv)