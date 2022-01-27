#!/usr/bin/env python3

import argparse, sys, re

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="FOLD", metavar="FOLD", default=70,type= int)
parser.add_argument("-i", help="Input file name", metavar="[filename]", required=True, type=str)

args = parser.parse_args()

fold = args.f
fileEMBL = open(args.i, "r")
filename = args.i.split(".",1)[0]


def file_checker_and_reader(inputfile):
    fileinlinelist = []
    for line in inputfile:
        line = line.rstrip()
        fileinlinelist.append(line)

    for line in fileinlinelist:
        if re.match("^ID",line) != None:
            return f"EMBL",fileinlinelist
        elif re.match("^#MEGA", line) != None:
            return f"MEGA",fileinlinelist
        elif re.match("^LOCUS", line) != None:
            return f"GenBank",fileinlinelist
        elif re.match("^@", line) != None:
            if re.match("^@SQ", line) != None or re.match("^@PQ", line) != None:
                return  f"SAM",fileinlinelist
            else:
                return f"FASTQ",fileinlinelist
        elif re.match("^##fileformat",line) != None:
            return f"VCF",fileinlinelist


fileformat, readablefilelist = file_checker_and_reader(fileEMBL)

class FASTA:
    def __init__(self):
        self.seq = []
        self.desc = []

def print_ready_seq(seq):
    outseq = ""

    for index in range(int(len(seq)/fold)+1):
        if index >= len(seq)/fold:
            outseq + "\n" + seq[index*fold:]
        else:
            outseq = outseq + "\n" + seq[index*fold:fold*(index+1)]
    outseq = outseq.lstrip()
    return outseq

output = FASTA()

def conversionToFASTA(fileformat,readablefilelist):

    #EMBL Operations
    if fileformat == "EMBL":
        ctr = 1
        seqctr = 0
        for line in readablefilelist:
            if re.match("^AC",line) != None:
                line = line.rstrip()
                line = line.split(" ",1)
                line[1] = line[1].lstrip()
                output.desc.append(line[1] + "|")
                output.seq.append("")
                seqctr += 1
            elif re.match("^DE",line) != None:
                line = line.split(" ",1)
                line[1] = line[1].lstrip()
                output.desc[seqctr-1] = output.desc[seqctr-1] + line[1]
            elif re.match("^SQ",line) != None:
                startseq = ctr
                for index in range(startseq, len(readablefilelist) - 1):
                    seq = readablefilelist[index]
                    if re.match("//", seq) != None:
                        break
                    seq = seq.strip()
                    seq = seq.split("  ", 1)
                    seq[0] = seq[0].upper()
                    output.seq[seqctr - 1] = output.seq[seqctr - 1] + seq[0]
                    output.seq[seqctr - 1] = output.seq[seqctr - 1].replace(" ", "")
            ctr += 1


    # GenBank operations
    elif fileformat == "GenBank" :
        ctr=1
        seqctr=0
        for line in readablefilelist:
            if re.match("^LOCUS",line) != None:
                line = line.replace(" ", "\t",1).split("\t")
                line[1] = line[1].strip()
                line[1] = line[1].split(" ",1)
                output.desc.append(line[1][0] + " ")
                output.seq.append("")
                seqctr += 1
            elif re.match("^DEF", line) != None:
                line = line.split(" ",1)
                line[1] = line[1].strip()
                output.desc[seqctr-1] = output.desc[seqctr-1] + line[1]
            elif re.match("^ORIGIN", line) != None:
                seqstart = ctr
                for index in range(seqstart, len(readablefilelist) - 1):
                    seq = readablefilelist[index]
                    if re.match("//", seq) != None:
                        break
                    seq = seq.strip()
                    seq = seq.split(" ", 1)
                    seq[1] = seq[1].upper()
                    output.seq[seqctr - 1] = output.seq[seqctr - 1] + seq[1]
                    output.seq[seqctr - 1] = output.seq[seqctr - 1].replace(" ", "")
            ctr += 1

    # MEGA operations
    if fileformat == "MEGA":
        ctr = 1
        seqctr = 0
        for line in readablefilelist:
            if re.match("^#MEGA", line) != None or re.match("^#mega", line) != None:
                continue
            elif re.match("^TITLE", line) != None:
                continue
            elif re.match("^#", line) != None:
                line = line.replace("#","",1)
                output.desc.append(line)
                output.seq.append("")
                seqctr += 1
            elif re.match("[ACGTNacgtn]", line) != None:
                seq = re.match("[ACGTNacgtn]+",line).string
                output.seq[seqctr-1] = output.seq[seqctr-1] + seq

    # FASTQ Operations
    if fileformat == "FASTQ":
        seqctr = 0
        for index in range(0, len(readablefilelist)):
            line = readablefilelist[index]
            if re.match("^@", line) != None:
                line = line.strip()
                line = line.replace("@","",1)
                output.desc.append(line)
                output.seq.append("")
                seqctr += 1
                continue
            if index > 0:
                prevline = readablefilelist[index - 1]
                if re.match("^\+", line):
                    continue
                elif re.match("^\+", prevline):
                    continue
            if re.match("^\n", line) == None or re.match("^@", line) == None or re.match("^\+", line) == None:
                line = line.strip()
                output.seq[seqctr - 1] = output.seq[seqctr - 1] + line

    # SAM Operations
    if fileformat == "SAM":
        seqctr = 0
        for index in range(0,len(readablefilelist)):
            line = readablefilelist[index]
            if re.match("^@", line) != None:
                continue
            else:
                line = line.split("\t")
                output.desc.append(line[0])
                output.seq.append(line[9])

    # VCF Operations
    if fileformat == "VCF":
        seqctr = 0
        refDict = {}
        output.seq = {}
        storeindex = 0
        for index in range(0, len(readablefilelist)):
            line = readablefilelist[index]
            if re.match("^##", line) != None:
                continue
            elif re.match("^#CHROM", line) != None:
                line = line.split()
                ctr = 0
                for item in line:
                    refDict[item] = []
                    ctr += 1
                    if ctr >= 10:
                        output.desc.append(item)
            else:
                line = line.split("\t")
                ctr = 0
                if not line[0] in output.desc:
                    output.desc.insert(seqctr,line[0])
                    storeindex =index
                    seqctr += 1
                for item in refDict:
                    refDict[item].append(line[ctr])
                    ctr += 1

        for desc in output.desc:
            output.seq[desc] = ""
            if desc not in refDict.keys():
                for seq in refDict['REF']:
                    if desc == refDict['#CHROM'][refDict['REF'].index(seq)]:
                        output.seq[desc] = output.seq[desc] + seq
            else:
                for index in range(0, len(refDict[desc])):
                    value = refDict[desc][index]
                    value = value.split(":",1)
                    if value[0] == '0':
                        output.seq[desc] = output.seq[desc] + refDict['REF'][index]
                    else:
                        # seq = refDict['ALT'][index]
                        output.seq[desc] = output.seq[desc] + refDict['ALT'][index].split(",")[int(value[0])-1]
        output.seq = list(output.seq.values())

conversionToFASTA(fileformat,readablefilelist)

def FASTAprinter(desc, seq, file):
    for index in range(0, len(desc)):
        file.write(f">{desc[index]}\n{print_ready_seq(seq[index])}\n")

if re.fullmatch("[ACGTNacgtn]", output.seq[0]) is not None:
    residue = ".faa"
else:
    residue = ".fna"

outputfile = open(filename+residue,"w")
FASTAprinter(output.desc, output.seq, outputfile)