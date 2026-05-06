#!/usr/bin/python3
from collections import OrderedDict
import sys
import getopt
import re

USAGE = "\n\nThe script summarizes statistic of chemotaxis systems created by the process_MiST_ST_counts.py script.\n" + \
    "For each taxonomy level (species, genus, family, order, class, phylum, kingdom) it calculates the average number of each chemotaxis system component per genome.\n" + \
    "python 	" + sys.argv[0] + '''
    -h || --help               - help
    -c || --cfile              - input file with cemotaxis data  
    -m || --mfile              - input file MCP chemoreceptors
    -s || --sfile              - input file (GTDB taxonomy metadata file)
    -t || --taxlevel           - taxonomy level for summarization. One of: species, genus, family, order, class, taxlevel, kingdom, or acorss. 'across' means across all phyla.

    Input files are tab delimited and have the following format (ex., -a file, CheA statistics file):
    genome	genome_accession	protein_type	total	chew	checx	other	MAC1	MAC2	F1	F2	F3	F4	F5	F6	F7	F8	F9	F10	F11	F12	F13	F14	F15	F16	F17	ACF	Tfp	total_norm_by_protein_count	chew_norm_by_protein_count	checx_norm_by_protein_count	other_norm_by_protein_count	MAC1_norm_by_protein_count	MAC2_norm_by_protein_count	F1_norm_by_protein_count	F2_norm_by_protein_count	F3_norm_by_protein_count	F4_norm_by_protein_count	F5_norm_by_protein_count	F6_norm_by_protein_count	F7_norm_by_protein_count	F8_norm_by_protein_count	F9_norm_by_protein_count	F10_norm_by_protein_count	F11_norm_by_protein_count	F12_norm_by_protein_count	F13_norm_by_protein_count	F14_norm_by_protein_count	F15_norm_by_protein_count	F16_norm_by_protein_count	F17_norm_by_protein_count	ACF_norm_by_protein_count	Tfp_norm_by_protein_count
    GCA_001872605.1	GCA_001872605	chea	1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00021272069772388852	0	0	0	0	0	0.00021272069772388852	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    GCA_001872605.1	GCA_001872605	chew	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00021272069772388852	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    GCA_001872605.1	GCA_001872605	checx	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.00042544139544777704	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    '''
CHEM_TEMPLATE_TO_COUNTS = OrderedDict([("total", 0), ("chew", 0), ("checx", 0), ("other", 0), ("MAC1", 0), ("MAC2", 0), ("F1", 0), ("F2", 0), ("F3", 0), ("F4", 0), ("F5", 0), ("F6", 0), ("F7", 0), ("F8", 0), ("F9", 0), ("F10", 0), ("F11", 0), ("F12", 0), ("F13", 0), ("F14", 0), ("F15", 0), ("F16", 0), ("F17", 0), ("ACF", 0), ("Tfp", 0), \
								("total_norm_by_protein_count", 0), ("chew_norm_by_protein_count", 0), ("checx_norm_by_protein_count", 0), ("other_norm_by_protein_count", 0), ("MAC1_norm_by_protein_count", 0), ("MAC2_norm_by_protein_count", 0), ("F1_norm_by_protein_count", 0), ("F2_norm_by_protein_count", 0), \
								("F3_norm_by_protein_count", 0), ("F4_norm_by_protein_count", 0), ("F5_norm_by_protein_count", 0), ("F6_norm_by_protein_count", 0), ("F7_norm_by_protein_count", 0), ("F8_norm_by_protein_count", 0), ("F9_norm_by_protein_count", 0), \
								("F10_norm_by_protein_count", 0), ("F11_norm_by_protein_count", 0), ("F12_norm_by_protein_count", 0), ("F13_norm_by_protein_count", 0), ("F14_norm_by_protein_count", 0), ("F15_norm_by_protein_count", 0), \
								("F16_norm_by_protein_count", 0), ("F17_norm_by_protein_count", 0), ("ACF_norm_by_protein_count", 0), ("Tfp_norm_by_protein_count", 0), ("recordNumber", 0)])

CHEA_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEW_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHECX_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEMOTHER_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHED_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEZ_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEV_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHEB_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()
CHER_TO_COUNTS = CHEM_TEMPLATE_TO_COUNTS.copy()

MCP_TO_COUNTS = OrderedDict([("total", 0), ("64H", 0), ("58H", 0), ("52H", 0), ("48H", 0), ("44H", 0), ("42H", 0), ("40H", 0), ("38H", 0), ("36H", 0), ("34H", 0), ("28H", 0), ("24H", 0), \
                             	("total_norm_by_protein_count", 0), ("64H_norm_by_protein_count", 0), ("58H_norm_by_protein_count", 0), ("52H_norm_by_protein_count", 0), \
								("48H_norm_by_protein_count", 0), ("44H_norm_by_protein_count", 0), ("42H_norm_by_protein_count", 0), ("40H_norm_by_protein_count", 0), \
								("38H_norm_by_protein_count", 0), ("36H_norm_by_protein_count", 0), ("34H_norm_by_protein_count", 0), ("28H_norm_by_protein_count", 0), ("24H_norm_by_protein_count", 0), ("recordNumber", 0)]) 

CHEM_COMPONENT_TO_DATASTRUCTURE = dict([("chea", CHEA_TO_COUNTS), ("chew", CHEW_TO_COUNTS), ("checx", CHECX_TO_COUNTS), ("other", CHEMOTHER_TO_COUNTS), ("chev", CHEV_TO_COUNTS), ("cheb", CHEB_TO_COUNTS), ("cher", CHER_TO_COUNTS), ("ched", CHED_TO_COUNTS), ("chez", CHEZ_TO_COUNTS)]) 
MCP_COMPONENT_TO_DATASTRUCTURE = dict([("mcp", MCP_TO_COUNTS)]) 

GTDB_FILE = None
CHEM_FILE = None
MCP_FILE = None

GENOME_TO_CHEA_NUMBER = {}

TAXLEVEL = "phylum"
TAXONOMY_TO_LEVEL = {"species": 7, "genus": 6, "family": 5, "order": 4, "class": 3, "phylum": 2, "kingdom": 1}
GENOME_TO_TAXONOMY = {}
ADDITIONAL_HEADERS = "\t".join(["gtdb_taxonomy_string", "gtdb_taxonomy_last", "gtdb_taxonomy_rank"])

def initialize(argv):
    global INPUT_FILE_TO_DATA, CHEM_FILE, MCP_FILE, GTDB_FILE, TAXLEVEL
    try:
        opts, args = getopt.getopt(argv[1:],"hi:c:m:s:t:",["help", "cfile=", "mfile", "sfile=", "taxlevel="])
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
            elif opt in ("-c", "--afile"):
                CHEM_FILE = str(arg).strip()
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

def processSTstatistics(fileToProcess, componen_to_datastructure, mcp=False):
    taxlevel_to_data = {}
    with open (fileToProcess, "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            if lineNumber > 0:
                line = line.strip().split("\t")
                taxlevel = GENOME_TO_TAXONOMY[line[0]]
                if lineNumber == 1 or (lineNumber > 1 and taxlevel not in taxlevel_to_data):
                    taxlevel_to_data[taxlevel] = componen_to_datastructure.copy()

                if not mcp:
                    protein_type = line[1]
                    taxlevel_to_data[taxlevel][protein_type]["recordNumber"] += 1
                    # I use idx+shift becasue the first 1 or 2 columns in the input file are extra conpared to the CHEM_TEMPLATE_TO_COUNTS & MCP_TO_COUNTS fields.
                    # Particularly, in the input file of chemo systems (CHEM_FILE) there are two additional fields: genome identifier, and potein type
                    # while in the input file of MCP chemoreceptors (MCP_FIKE) there is one additional field: genome identifier
                    # while elements of CHEM_TEMPLATE_TO_COUNTS/MCP_TO_COUNTS look like: "protein_type", "total", "F1, F2", etc for chemo systems, or "total", "64H", etc. for mcp chemoreceptors 
                    # Therefore we add 'shift' value to shift to correspinding positions from the beginning of the line and land on the next column. Look at the input file format in USAGE
                    shift = 2
                    for idx, component in enumerate(taxlevel_to_data[taxlevel][protein_type]):
                        if component != "recordNumber": 
                            taxlevel_to_data[taxlevel][protein_type][component] += float(line[idx+shift])
                else:
                    taxlevel_to_data[taxlevel]["mcp"]["recordNumber"] += 1
                    shift = 1
                    for idx, component in enumerate(taxlevel_to_data[taxlevel]["mcp"]):
                        if component != "recordNumber": 
                            taxlevel_to_data[taxlevel]["mcp"][component] += float(line[idx+shift])              

    return taxlevel_to_data

def finalizeDataAndSave(taxlevel_to_data, inputFile, data, mcp=False):
    protein_type_header = "protein_type" if not mcp else ""
    with open (inputFile.split(".")[0]+"_"+TAXLEVEL, "a") as output_file:
         # converting to list and write all the headers except for the last one, which is 'record_number'
        output_file.write(ADDITIONAL_HEADERS + "\t" + protein_type_header + "\t" + "\t".join(list(CHEM_TEMPLATE_TO_COUNTS.keys())[:-1]) + "\n")

        for taxlevel, data in taxlevel_to_data.items():
            for protein_type, chem_mcp_info_data in data.items():
                for component in chem_mcp_info_data:
                    chem_mcp_info_data[component] = chem_mcp_info_data[component]/chem_mcp_info_data["recordNumber"]
                output_file.write("\t".join([taxlevel, taxlevel.split(";")[-1], TAXLEVEL, protein_type if not mcp else ""])  + "\t" + "\t".join(list(map(str, chem_mcp_info_data.values()))[:-1]) + "\n")

def main(argv):
    initialize(argv)
    # Process chemotaxis components
    taxlevel_to_data = processSTstatistics(CHEM_FILE, CHEM_COMPONENT_TO_DATASTRUCTURE, False)
    finalizeDataAndSave(taxlevel_to_data, CHEM_FILE, CHEM_COMPONENT_TO_DATASTRUCTURE, False)
    
    # Process MCP chemoreceptors
    taxlevel_to_data = processSTstatistics(MCP_FILE, MCP_COMPONENT_TO_DATASTRUCTURE, True)
    finalizeDataAndSave(taxlevel_to_data, MCP_FILE, MCP_COMPONENT_TO_DATASTRUCTURE, True)

        

main(sys.argv)