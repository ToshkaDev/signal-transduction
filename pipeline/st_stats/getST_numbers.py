#!/usr/bin/python3
import sys
from collections import OrderedDict

#PHYLUM_TO_DATA = {"pylum1": {"total": 0, ...}}
PHYLUM_TO_DATA = {}
DATA = OrderedDict([("total", 0), ("ocpTotal", 0), ("tcpTotal", 0), ("tcpHKTotal", 0), ("tcpRRTotal", 0), ("tcpHK", 0), ("tcpHHK", 0), ("tcpRR", 0), ("tcpHRR", 0), \
                    ("tcpTotal/ocpTotal", 0), ("tcpHKTotal/tcpRRTotal", 0), ("tcpHK/tcpRR", 0), ("tcpHHK/tcpHRR", 0), \
                    ("tcpOther", 0), ("chemSys", 0), ("ecf", 0), ("other", 0), ("recordNumber", 0)])

#GenomeId(0)    Total(1)   OCP_total(2)    TCP_total(3) TCP_HK(4)   TCP_HHK(5)  TCP_RR(6)   TCP_HRR(7)  TCP_Other(8)    ChemSys(9) ECF(10) Other(11)   Taxonomy-Phylum(12)
def processSTstatistics():
    with open (sys.argv[1], "r") as inputFile:
        for lineNumber, line in enumerate(inputFile):
            if lineNumber > 0:
                line = line.strip().split("\t")
                if len(line) >= 13:
                    phylum = line[12]
                else:
                    phylum = "Across_Phyla"
                if lineNumber == 1:
                    PHYLUM_TO_DATA[phylum] = DATA.copy()
                elif phylum not in PHYLUM_TO_DATA and lineNumber > 1:
                    PHYLUM_TO_DATA[phylum] = DATA.copy()
                    
                PHYLUM_TO_DATA[phylum]["recordNumber"] += 1
                PHYLUM_TO_DATA[phylum]["total"] += int(line[1])
                PHYLUM_TO_DATA[phylum]["ocpTotal"] += int(line[2])
                PHYLUM_TO_DATA[phylum]["tcpTotal"] += int(line[3])
                PHYLUM_TO_DATA[phylum]["tcpHK"] += int(line[4])
                PHYLUM_TO_DATA[phylum]["tcpHHK"] += int(line[5])
                PHYLUM_TO_DATA[phylum]["tcpRR"] += int(line[6])
                PHYLUM_TO_DATA[phylum]["tcpHRR"] += int(line[7])
                PHYLUM_TO_DATA[phylum]["tcpOther"] += int(line[8])
                PHYLUM_TO_DATA[phylum]["chemSys"] += int(line[9])
                PHYLUM_TO_DATA[phylum]["ecf"] += int(line[10])
                PHYLUM_TO_DATA[phylum]["other"] += int(line[11])
                PHYLUM_TO_DATA[phylum]["tcpHKTotal"] += (int(line[4]) + int(line[5]))
                PHYLUM_TO_DATA[phylum]["tcpRRTotal"] += (int(line[6]) + int(line[7]))
                if float(line[2]):
                    PHYLUM_TO_DATA[phylum]["tcpTotal/ocpTotal"] += int(line[3])/float(line[2])
                if float(line[6]):
                    PHYLUM_TO_DATA[phylum]["tcpHK/tcpRR"] += int(line[4])/float(line[6])
                if float(line[7]):
                    PHYLUM_TO_DATA[phylum]["tcpHHK/tcpHRR"] += int(line[5])/float(line[7])
                if float(PHYLUM_TO_DATA[phylum]["tcpRRTotal"]):
                    PHYLUM_TO_DATA[phylum]["tcpHKTotal/tcpRRTotal"] += PHYLUM_TO_DATA[phylum]["tcpHKTotal"]/float(PHYLUM_TO_DATA[phylum]["tcpRRTotal"])
                    
def finalizeDataAndPrint():
    for phylum, data in PHYLUM_TO_DATA.items():
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
        print(phylum + "\t" + "\t".join(map(str, dataValues)))

def roundToFirstDecim(value):
    return round(value, 1)
	
def main():
    print ("Phylum" + "\t" + "\t".join(DATA.keys()))
    processSTstatistics()
    finalizeDataAndPrint()

main()