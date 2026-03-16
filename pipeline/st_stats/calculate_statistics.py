#!/usr/bin/python3
import sys
import getopt
import re
from collections import OrderedDict

USAGE = "\n\nThe script summarizes statistic created by the process_MiST_ST_counts.py script.\n\n" + \
    "python 	" + sys.argv[0] + '''
    -h || --help               - help
    -i || --ifile              - input file (statistics file). Content example:
        genomeID	 Total	 OCP_total	 TCP_total	 TCP_HK	 TCP_HHK	 TCP_RR	 TCP_HRR	 TCP_Other	 ChemSys	 ECF	 Other	 Taxonomy	 
        GCA_001872605.1	201	72	100	30	4	58	6	2	1	8	14	d__Bacteria;p__Armatimonadota;c__Abditibacteria;o__CG2-30-59-28;f__CG2-30-59-28;g__CG2-30-59-28;s__CG2-30-59-28 sp001872605
        GCA_001873295.1	145	43	77	25	9	37	6	0	2	3	11	d__Bacteria;p__CG2-30-70-394;c__CG2-30-70-394;o__CG2-30-70-394;f__CG2-30-70-394;g__CG2-30-70-394;s__CG2-30-70-394 sp001873295
    -s || --sfile              - input file 2 (GTDB taxonomy metadata file)
    -t || --taxlevel           - taxonomy level for summarization. One of: species, genus, family, order, class, taxlevel, kingdom, or acorss. 'across' meas across all phyla
    -o || --ofile              - output file
    '''

INPUT_FILE = None
INPUT_FILE2 = None
OUTPUT_FILE = "st_stats_summary.txt"
TAXLEVEL = "phylum"

TAXONOMY_TO_LEVEL = {"species": 7, "genus": 6, "family": 5, "order": 4, "class": 3, "phylum": 2, "kingdom": 1}
#TAXONOMY_LEVEL_TO_DATA = {"pylum1": {"total": 0, ...}}
TAXONOMY_LEVEL_TO_DATA = {}
DATA = OrderedDict([("total", 0), ("ocp_total", 0), ("tcp_total", 0), ("tcp_hk_total", 0), ("tcp_rr_total", 0), ("tcp_hk", 0), ("tcp_hhk", 0), ("tcp_rr", 0), ("tcp_hrr", 0), \
                    ("tcp_total_by_ocp_total", 0), ("tcp_hk_total_by_tcp_rr_total", 0), ("tcp_hk_by_tcp_rr", 0), ("tcp_hhk_by_tcp_hrr", 0), \
                    ("tcp_other", 0), ("chem_sys", 0), ("ecf", 0), ("other", 0), ("record_number", 0)])
GENOME_TO_TAXONOMY = {}
ADDITIONAL_HEADERS = "\t".join(["gtdb_taxonomy_string", "gtdb_taxonomy_last", "gtdb_taxonomy_rank"])

def initialize(argv):
    global INPUT_FILE, INPUT_FILE2, OUTPUT_FILE, TAXLEVEL
    try:
        opts, args = getopt.getopt(argv[1:],"hi:s:o:t:",["help", "ifile=", "sfile=", "ofile=", "taxlevel="])
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
            elif opt in ("-s", "--sfile"):
                INPUT_FILE2 = str(arg).strip()
            elif opt in ("-o", "--ofile"):
                OUTPUT_FILE = str(arg).strip()
            elif opt in ("-t", "--taxlevel"):
                TAXLEVEL = str(arg).strip().lower()
                OUTPUT_FILE += "_" + TAXLEVEL
    except Exception as e:
        print(("===========ERROR==========\n " + str(e) + USAGE))
        sys.exit(2)
            
def processSTstatistics():
    # regex to remove prefix .__ (like d__, p__, etc. from the GTDB taxonomomy string)
    regex=r'.__'
    with open(INPUT_FILE2, "r") as iFile2:
        for line in iFile2:
            #$0 - genome version, $1 - genome accession, $2 - genome size, $3 - protein counts, $4 - GTDB taxonomy, $5 - NCBI taxonomy
            record = line.strip().split("\t")
            genome_version = record[0]
            taxonomy = ";".join(record[4].split(";")[:TAXONOMY_TO_LEVEL[TAXLEVEL]])
            taxonomy = re.sub(regex, "", taxonomy)
            GENOME_TO_TAXONOMY[genome_version] = taxonomy

    # Main logic
    with open (INPUT_FILE, "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            if lineNumber > 0:
                line = line.strip().split("\t")
                taxlevel = GENOME_TO_TAXONOMY[line[0]]

                if lineNumber == 1 or (lineNumber > 1 and taxonomy not in TAXONOMY_LEVEL_TO_DATA):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel] = DATA.copy()

                TAXONOMY_LEVEL_TO_DATA[taxlevel]["record_number"] += 1
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["total"] += int(line[2])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["ocp_total"] += int(line[3])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_total"] += int(line[4])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk"] += int(line[5])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hhk"] += int(line[6])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr"] += int(line[7])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hrr"] += int(line[8])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_other"] += int(line[9])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["chem_sys"] += int(line[10])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["ecf"] += int(line[11])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["other"] += int(line[12])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total"] += (int(line[5]) + int(line[6]))
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total"] += (int(line[7]) + int(line[8]))
                if float(line[3]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_total_by_ocp_total"] += int(line[4])/float(line[3])
                if float(line[7]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_by_tcp_rr"] += int(line[5])/float(line[7])
                if float(line[8]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hhk_by_tcp_hrr"] += int(line[6])/float(line[8])
                if float(TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total"]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total_by_tcp_rr_total"] += TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total"]/float(TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total"])
                    
def finalizeDataAndPrint():
    #Write headers first
    with open(OUTPUT_FILE, 'a') as output_file:
            output_file.write(ADDITIONAL_HEADERS + "\t" + "\t".join(DATA.keys()) + "\n")
            
    for taxlevel, data in TAXONOMY_LEVEL_TO_DATA.items():
        data["total"] = data["total"]/data["record_number"]
        data["ocp_total"] = data["ocp_total"]/data["record_number"]
        data["tcp_total"] = data["tcp_total"]/data["record_number"]
        data["tcp_hk"] = data["tcp_hk"]/data["record_number"]
        data["tcp_hhk"] = data["tcp_hhk"]/data["record_number"]
        data["tcp_rr"] = data["tcp_rr"]/data["record_number"]
        data["tcp_hrr"] = data["tcp_hrr"]/data["record_number"]
        data["tcp_other"] = data["tcp_other"]/data["record_number"]
        data["chem_sys"] = data["chem_sys"]/data["record_number"]
        data["ecf"] = data["ecf"]/data["record_number"]
        data["other"] = data["other"]/data["record_number"]
        data["tcp_hk_total"] = data["tcp_hk_total"]/data["record_number"]
        data["tcp_rr_total"] = data["tcp_rr_total"]/data["record_number"]      
        data["tcp_total_by_ocp_total"] = data["tcp_total_by_ocp_total"]/data["record_number"]
        data["tcp_hk_by_tcp_rr"] = data["tcp_hk_by_tcp_rr"]/data["record_number"]
        data["tcp_hhk_by_tcp_hrr"] = data["tcp_hhk_by_tcp_hrr"]/data["record_number"]
        data["tcp_hk_total_by_tcp_rr_total"] = data["tcp_hk_total_by_tcp_rr_total"]/data["record_number"]
        dataValues = map(roundToFirstDecim, data.values())
        #And now write data
        with open(OUTPUT_FILE, 'a') as output_file:
            output_file.write("\t".join([taxlevel, taxlevel.split(";")[-1], TAXLEVEL]) + "\t" + "\t".join(map(str, dataValues)) + "\n")

def roundToFirstDecim(value):
    return round(value, 1)
    
def main(argv):
    initialize(argv)
    processSTstatistics()
    finalizeDataAndPrint()

main(sys.argv)