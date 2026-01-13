#!/usr/bin/python3
import sys
import getopt
from collections import OrderedDict

USAGE = "\n\nThe script summarizes statistic created by the process_MiST_ST_counts.py script.\n\n" + \
    "python 	" + sys.argv[0] + '''
    -h || --help               - help
    -i || --ifile              - input file (statistics file). Content example:
        genomeID	 Total	 OCP_total	 TCP_total	 TCP_HK	 TCP_HHK	 TCP_RR	 TCP_HRR	 TCP_Other	 ChemSys	 ECF	 Other	 Taxonomy	 
        GCA_001872605.1	201	72	100	30	4	58	6	2	1	8	14	d__Bacteria;p__Armatimonadota;c__Abditibacteria;o__CG2-30-59-28;f__CG2-30-59-28;g__CG2-30-59-28;s__CG2-30-59-28 sp001872605
        GCA_001873295.1	145	43	77	25	9	37	6	0	2	3	11	d__Bacteria;p__CG2-30-70-394;c__CG2-30-70-394;o__CG2-30-70-394;f__CG2-30-70-394;g__CG2-30-70-394;s__CG2-30-70-394 sp001873295
    -t || --taxlevel           - taxonomy level for summarization. One of: species, genus, family, order, class, taxlevel, kingdom, or acorss. 'across' meas across all phyla
    -o || --ofile              - output file
    '''

INPUT_FILE = None
OUTPUT_FILE = "st_stats_summary.txt"
TAXLEVEL = "phylum"

TAXONOMY_TO_LEVEL = {"species": 7, "genus": 6, "family": 5, "order": 4, "class": 3, "phylum": 2, "kingdom": 1}
#TAXONOMY_LEVEL_TO_DATA = {"pylum1": {"total": 0, ...}}
TAXONOMY_LEVEL_TO_DATA = {}
DATA = OrderedDict([("total", 0), ("ocpTotal", 0), ("tcpTotal", 0), ("tcpHKTotal", 0), ("tcpRRTotal", 0), ("tcpHK", 0), ("tcpHHK", 0), ("tcpRR", 0), ("tcpHRR", 0), \
                    ("tcpTotal/ocpTotal", 0), ("tcpHKTotal/tcpRRTotal", 0), ("tcpHK/tcpRR", 0), ("tcpHHK/tcpHRR", 0), \
                    ("tcpOther", 0), ("chemSys", 0), ("ecf", 0), ("other", 0), ("recordNumber", 0)])


def initialize(argv):
    global INPUT_FILE, OUTPUT_FILE, TAXLEVEL
    try:
        opts, args = getopt.getopt(argv[1:],"hi:o:t:",["help", "ifile=", "ofile=", "taxlevel="])
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
            elif opt in ("-o", "--ofile"):
                OUTPUT_FILE = str(arg).strip()
            elif opt in ("-t", "--taxlevel"):
                TAXLEVEL = str(arg).strip().lower()
                OUTPUT_FILE += "_" + TAXLEVEL
    except Exception as e:
        print(("===========ERROR==========\n " + str(e) + USAGE))
        sys.exit(2)
            
def processSTstatistics():
    with open (INPUT_FILE, "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            if lineNumber > 0:
                line = line.strip().split("\t")
                if TAXLEVEL in TAXONOMY_TO_LEVEL:
                    taxlevel = ";".join(line[-1].split(";")[:TAXONOMY_TO_LEVEL[TAXLEVEL]])
                else:
                    taxlevel = "Across_Phyla"
                if lineNumber == 1 or (lineNumber > 1 and taxlevel not in TAXONOMY_LEVEL_TO_DATA):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel] = DATA.copy()
                    
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["recordNumber"] += 1
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["total"] += int(line[1])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["ocpTotal"] += int(line[2])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpTotal"] += int(line[3])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHK"] += int(line[4])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHHK"] += int(line[5])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpRR"] += int(line[6])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHRR"] += int(line[7])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpOther"] += int(line[8])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["chemSys"] += int(line[9])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["ecf"] += int(line[10])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["other"] += int(line[11])
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHKTotal"] += (int(line[4]) + int(line[5]))
                TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpRRTotal"] += (int(line[6]) + int(line[7]))
                if float(line[2]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpTotal/ocpTotal"] += int(line[3])/float(line[2])
                if float(line[6]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHK/tcpRR"] += int(line[4])/float(line[6])
                if float(line[7]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHHK/tcpHRR"] += int(line[5])/float(line[7])
                if float(TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpRRTotal"]):
                    TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHKTotal/tcpRRTotal"] += TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpHKTotal"]/float(TAXONOMY_LEVEL_TO_DATA[taxlevel]["tcpRRTotal"])
                    
def finalizeDataAndPrint():
    #Write headers first
    with open(OUTPUT_FILE, 'a') as output_file:
            output_file.write("taxLevel" + "\t" + "\t".join(DATA.keys()) + "\n")
            
    for taxlevel, data in TAXONOMY_LEVEL_TO_DATA.items():
        data["total"] = data["total"]/data["recordNumber"]
        data["ocpTotal"] = data["ocpTotal"]/data["recordNumber"]
        data["tcpTotal"] = data["tcpTotal"]/data["recordNumber"]
        data["tcpHK"] = data["tcpHK"]/data["recordNumber"]
        data["tcpHHK"] = data["tcpHHK"]/data["recordNumber"]
        data["tcpRR"] = data["tcpRR"]/data["recordNumber"]
        data["tcpHRR"] = data["tcpHRR"]/data["recordNumber"]
        data["tcpOther"] = data["tcpOther"]/data["recordNumber"]
        data["chemSys"] = data["chemSys"]/data["recordNumber"]
        data["ecf"] = data["ecf"]/data["recordNumber"]
        data["other"] = data["other"]/data["recordNumber"]
        data["tcpHKTotal"] = data["tcpHKTotal"]/data["recordNumber"]
        data["tcpRRTotal"] = data["tcpRRTotal"]/data["recordNumber"]      
        data["tcpTotal/ocpTotal"] = data["tcpTotal/ocpTotal"]/data["recordNumber"]
        data["tcpHK/tcpRR"] = data["tcpHK/tcpRR"]/data["recordNumber"]
        data["tcpHHK/tcpHRR"] = data["tcpHHK/tcpHRR"]/data["recordNumber"]
        data["tcpHKTotal/tcpRRTotal"] = data["tcpHKTotal/tcpRRTotal"]/data["recordNumber"]
        dataValues = map(roundToFirstDecim, data.values())
        #And now write data
        with open(OUTPUT_FILE, 'a') as output_file:
            output_file.write(taxlevel + "\t" + "\t".join(map(str, dataValues)) + "\n")

def roundToFirstDecim(value):
    return round(value, 1)
    
def main(argv):
    initialize(argv)
    processSTstatistics()
    finalizeDataAndPrint()

main(sys.argv)