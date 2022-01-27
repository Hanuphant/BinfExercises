#!/usr/bin/env python3

import sys

file1 = open(sys.argv[1],"r")
file1 = file1.read().strip()
seq1 = file1.split("\n")[1]
file2 = open(sys.argv[2],"r")
file2 = file2.read().strip()
seq2 = file2.split("\n")[1]

# print(seq1)?
# print(seq2)

seq1.upper()
seq2.upper()

# Sequence Checker
for nucl in list(seq1+seq2):
    if nucl=='A' or nucl=='G' or nucl=='T' or nucl=='C':
        continue
    else:
        print(f"Enter valid nucleotide sequence!!")
        exit(1)

matrix = [[ 0 for row in range(len(seq2)+2)] for column in range(len(seq1)+2)]
tracebackmatrix = [[ 0 for row in range(len(seq2)+2)] for column in range(len(seq1)+2)]

match = 1
mismatch = -1
gap = -1

# Initializing

# Label the sequences on the top and left side
for row in range(2, len(seq1)+2):
    matrix[row][0] = seq1[row-2]
    matrix[row][1] = matrix[row-1][1]+gap
for column in range(2, len(seq2)+2):
    matrix[0][column] = seq2[column-2]
    matrix[1][column] = matrix[1][column-1]+gap
def matrixprinter(mat):
    for row in mat:
        print(f"{row}")

def matcher(nucl1, nucl2):
    if nucl1 == nucl2:
        return match
    else:
        return mismatch

def scoreofcell(mat, row, column,mat2):
    score = max(mat[row-1][column-1]+matcher(mat[row][0],mat[0][column]), mat[row-1][column]+gap,mat[row][column-1]+gap)
    if score == mat[row-1][column-1]+matcher(mat[row][0],mat[0][column]):
        mat2[row][column] = "d"
    elif score == mat[row-1][column]+gap:
        mat2[row][column] = "u"
    elif score == mat[row][column-1]+gap:
        mat2[row][column] = "l"
    return score


for row in range(2, len(matrix)):
    for column in range(2, len(matrix[0])):
        matrix[row][column] = scoreofcell(matrix, row, column,tracebackmatrix)

# matrixprinter(tracebackmatrix)
# matrixprinter(matrix)

# Traceback
# print(f"\n")
traceback = []
def direction(mat, column, row):
    traceback.append(mat[row][column])
    if mat[row][column] == 'd':
        direction(mat,column-1,row-1)
    elif mat[row][column] == 'u':
        direction(mat,column,row-1)
    elif mat[row][column] == 'l':
        direction(mat, column-1,row)
    else:
        return 0

direction(tracebackmatrix, len(matrix[0])-1, len(matrix)-1)
traceback.pop()
traceback.reverse()
alignscore = 0
# matrixprinter(tracebackmatrix)
alignment = [[] for row in range(3)]
# matrixprinter(alignment)
# print(traceback)
ctr1 = 0
ctr2 = 0
for item in traceback:
        if ctr1<len(seq1) and ctr2<len(seq2):
            if item == 'd':
                if seq1[ctr1]==seq2[ctr2]:
                    alignment[1].append("|")
                    alignscore += match
                else:
                    alignment[1].append("*")
                    alignscore += mismatch
                alignment[0].append(seq1[ctr1])
                alignment[2].append(seq2[ctr2])
                ctr1+=1
                ctr2+=1
            elif item == 'l':
                alignment[1].append(" ")
                alignment[2].append(seq2[ctr2])
                alignment[0].append("-")
                alignscore += gap
                ctr2+=1
            elif item == 'u':
                alignment[1].append(" ")
                alignment[0].append(seq1[ctr1])
                alignment[2].append("-")
                alignscore += gap
                ctr1+=1

# matrixprinter(alignment)
# matrixprinter(matrix)
firstseq = ''
alignmentline = ''
secondseq = ''
print(f"{firstseq.join(alignment[0])}\n{alignmentline.join(alignment[1])}\n{secondseq.join(alignment[2])}")
print(f"Alignment Score : {alignscore}")