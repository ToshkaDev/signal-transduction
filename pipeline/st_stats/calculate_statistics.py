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
                    ("tcp_other", 0), ("chem_sys", 0), ("ecf", 0), ("other", 0), \
                    ("total_norm_by_avg_protein_count", 0), ("ocp_total_norm_by_avg_protein_count", 0), \
                    ("tcp_total_norm_by_avg_protein_count", 0), ("tcp_hk_total_norm_by_avg_protein_count", 0), \
                    ("tcp_rr_total_norm_by_avg_protein_count", 0), ("tcp_hk_total_norm_by_avg_protein_count", 0), \
                    ("tcp_hk_norm_by_avg_protein_count", 0), ("tcp_hhk_norm_by_avg_protein_count", 0), ("tcp_rr_norm_by_avg_protein_count", 0), \
                    ("tcp_hrr_norm_by_avg_protein_count", 0), ("tcp_total_by_ocp_total_norm_by_avg_protein_count", 0), \
                    ("tcp_hk_total_by_tcp_rr_total_norm_by_avg_protein_count", 0), ("tcp_hk_by_tcp_rr_norm_by_avg_protein_count", 0), \
                    ("tcp_hhk_by_tcp_hrr_norm_by_avg_protein_count", 0), ("tcp_other_norm_by_avg_protein_count", 0), \
                    ("chem_sys_norm_by_avg_protein_count", 0), ("ecf_norm_by_avg_protein_count", 0), ("other_norm_by_avg_protein_count", 0), ("record_number", 0)])

# To calculate statistic using helper tuple to avoid repeating code for each system. 
# Cannot use the DATA keys directly because the order of data in the input file is different from the order of keys in the DATA dictionary, 
# so we need to define a helper tuple that matches the order of data in the input file.
# System names in the helper tuple should match system names in the DATA dictionary, 
# but the order should match the order of data in the input file (starting from index 2, because the first two columns are genome_version and genome_accession).
HELPER_TUPLE = ("total", "ocp_total", "tcp_total", "tcp_hk", "tcp_hhk", "tcp_rr", "tcp_hrr", "tcp_other", "chem_sys", "ecf", "other", \
    "total_norm_by_avg_protein_count", "ocp_total_norm_by_avg_protein_count", "tcp_total_norm_by_avg_protein_count", \
"tcp_hk_norm_by_avg_protein_count", "tcp_hhk_norm_by_avg_protein_count", "tcp_rr_norm_by_avg_protein_count", \
"tcp_hrr_norm_by_avg_protein_count", "tcp_other_norm_by_avg_protein_count", "chem_sys_norm_by_avg_protein_count", \
"ecf_norm_by_avg_protein_count", "other_norm_by_avg_protein_count")

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
            
def processSTstatistics():
    # Main logic
    with open (INPUT_FILE, "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            if lineNumber > 0:
                line = line.strip().split("\t")
                taxlevel = GENOME_TO_TAXONOMY[line[0]]

                if lineNumber == 1 or (lineNumber > 1 and taxlevel not in TAXONOMY_LEVEL_TO_DATA):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel] = DATA.copy()

                TAXONOMY_LEVEL_TO_DATA[taxlevel]["record_number"] += 1
                # Using helper tuple to avoid repeating code for each system.
                for idx, system in enumerate(HELPER_TUPLE):
                    # +2 because the first two columns in the input file are genome_version and genome_accession, so the data starts from index 2
                    TAXONOMY_LEVEL_TO_DATA[taxlevel][system] += float(line[idx+2])
                
                # Compute combinations of systems and ratios
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total"] += (float(line[5]) + float(line[6]))
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total"] += (float(line[7]) + float(line[8]))
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total_norm_by_avg_protein_count"] += (float(line[16]) + float(line[17]))
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total_norm_by_avg_protein_count"] += (float(line[18]) + float(line[19]))

                if float(line[3]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_total_by_ocp_total"] += float(line[4])/float(line[3])
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_total_by_ocp_total_norm_by_avg_protein_count"] += float(line[15])/float(line[14])
                if float(line[7]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_by_tcp_rr"] += float(line[5])/float(line[7])
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_by_tcp_rr_norm_by_avg_protein_count"] += float(line[16])/float(line[18])
                if float(line[8]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hhk_by_tcp_hrr"] += float(line[6])/float(line[8])
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hhk_by_tcp_hrr_norm_by_avg_protein_count"] += float(line[17])/float(line[19])
                if float(TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total"]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total_by_tcp_rr_total"] += TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total"]/float(TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total"])
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total_by_tcp_rr_total_norm_by_avg_protein_count"] += TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_hk_total_norm_by_avg_protein_count"]/float(TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcp_rr_total_norm_by_avg_protein_count"])

def finalizeDataAndPrint():
    # Write headers first
    with open(OUTPUT_FILE, 'a') as output_file:
            # converting to list and write all the headers except for the last one, which is 'record_number'
            output_file.write(ADDITIONAL_HEADERS + "\t" + "\t".join(list(DATA.keys())[:-1]) + "\n")
    
    # Normalize data by record number
    for taxlevel, data in TAXONOMY_LEVEL_TO_DATA.items():
        for system in data.keys():
            data[system] = data[system]/data["record_number"]

        # And now write data
        with open(OUTPUT_FILE, 'a') as output_file:
            # Here converting the resulting map obejct to lost and writing all the values except for the last one (the number of records)
            output_file.write("\t".join([taxlevel, taxlevel.split(";")[-1], TAXLEVEL]) + "\t" + "\t".join(list(map(str, data.values()))[:-1]) + "\n")
   
def main(argv):
    initialize(argv)
    processSTstatistics()
    finalizeDataAndPrint()

main(sys.argv)