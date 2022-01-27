#!/usr/bin/env python3

import sys, re, subprocess, argparse
from multiprocessing import Pool

class Parallel_Ani:

    def __init__(self, files):
        self.files = files


    def unit_dna_diff(self,arguments):
        row, col, file1, file2 = arguments
        if row == col:
            align = 100
        else:
            output = str(row)+"_"+str(col)
            subprocess.check_output(['dnadiff','-p',output,file1,file2])
            align = self.fetch(output)
            subprocess.check_output(['rm',output+'.1coords'])
            subprocess.check_output(['rm',output+'.1delta'])
            subprocess.check_output(['rm',output+'.delta'])
            subprocess.check_output(['rm',output+'.mcoords'])
            subprocess.check_output(['rm',output+'.qdiff'])
            subprocess.check_output(['rm',output+'.mdelta'])
            subprocess.check_output(['rm',output+'.rdiff'])
            subprocess.check_output(['rm',output+'.report'])
            subprocess.check_output(['rm',output+'.snps'])
        return [row, col, align]

    @staticmethod
    def fetch(file):
        filer = open(file+".report","r")
        liner = filer.readlines()
        liner = liner[18]
        liner = liner.split()
        return float(liner[1])

if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument('-o', metavar="Output File", help="Output file", type=str, required=True)
    parser.add_argument('-t', metavar="Threads", help="Enter no of threads", type=int, required=False, default=5)
    parser.add_argument('fila', nargs = '+', type=str, help = "Enter files by spacing them" )
    args = parser.parse_args()
    
    filu = args.fila
    thread = args.t
    if thread > 5:
        thread = 5
        print(f"Warning. Threads passed were more than 5. For the safety of the program the threads are set default to 5.")

    filu = [''] + filu
    tabl = [[ 0  for column in range(len(filu)) ] for row in range(len(filu)) ]

    p = Parallel_Ani(filu)

    # Titles added
    for col in range(len(p.files)):
        tabl[0][col] = p.files[col]
        tabl[col][0] = p.files[col]

    # Preconditioning
    tasklist = []
    for col in range(1,len(p.files)):
        row = 1
        while row<=col:
            tasklist.append([row, col, tabl[0][col], tabl[row][0]])
            row += 1

    print(tasklist)

    pool = Pool(thread)
    result = list(pool.map(p.unit_dna_diff,tasklist))
    pool.close()
    pool.join()

    for element in result:
        tabl[element[0]][element[1]] = element[2]
        tabl[element[1]][element[0]] = element[2]

    outputfile = open(args.o,"w")
    for row in tabl:
        for element in row:
            outputfile.write(f"{element}\t")
        outputfile.write("\n")
