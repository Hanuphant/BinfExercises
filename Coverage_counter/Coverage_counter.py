#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help="Input File Name", metavar="[File Name]", required=True, type=str)

args = parser.parse_args()

BEDfile = open(args.i, "r")

readlist = [[],[],[]]
chromosome_index_dict = {}

for index,line in enumerate(BEDfile):
    line = line.strip()
    line = line.split("\t", 3)
    readlist[0].append(line[0])
    min_index, max_index = chromosome_index_dict.get(line[0], [10**8, -1])
    chromosome_index_dict[line[0]] = [min(index, min_index), max(index, max_index)]
    readlist[1].append(int(line[1]))
    readlist[2].append(int(line[2]))

def frequency_count(starts, ends):
    ans = []
    startpoints = [(start, "start") for start in starts]
    endpoints = [(end, "end") for end in ends]
    sorted_events = sorted(startpoints + endpoints)
    ctr = 1
    for index in range(1, len(sorted_events)):
        if ctr != 0 and sorted_events[index - 1][0] != sorted_events[index][0]:
            ans.append([sorted_events[index - 1][0], sorted_events[index][0], ctr])
        if sorted_events[index][1] == "start":
            ctr += 1
        else:
            ctr -= 1

    return ans

frequencies = {}
for chromosome_name, (start, end) in (chromosome_index_dict.items()):
    frequencies[chromosome_name] = frequency_count(readlist[1][start:end+1], readlist[2][start:end+1])
    for start_,end_,count in frequencies[chromosome_name]:
        print(f"{chromosome_name}\t{start_}\t{end_}\t{count}\t")
