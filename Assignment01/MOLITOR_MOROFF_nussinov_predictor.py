"""Dominik Molitor, Friederike Moroff, 30.04.2022"""
import sys
import numpy as np

def read_data_in():
    path = ""
    min_loop_len = 3
    score_gc = 1
    score_au = 1
    score_gu = 1
    for value in range(0, len(sys.argv)):
        if "-i" not in sys.argv:
            print("Error: Path is missing")
            quit()
        if sys.argv[value] == "-i":
            if not sys.argv[value+1].startswith("-"):
                path = sys.argv[value+1]
            else:
                print("Syntax Error: Incorrect syntax at -i")
                quit()
        elif sys.argv[value] == "--min-loop-length":
            try:
                min_loop_len = int(sys.argv[value+1])
            except:
                print("Syntax Exception: Invalid minimum loop length\nThe standard length of 3 will be used")
        elif sys.argv[value] == "--score-GC":
            try:
                score_gc = float(sys.argv[value+1])
            except:
                print("Syntax Exception: Invalid score value")
        elif sys.argv[value] == "--score-AU":
            try:
                score_au = float(sys.argv[value+1])
            except:
                print("Syntax Exception: Invalid score value")
        elif sys.argv[value] == "--score-GU":
            try:
                score_gu = float(sys.argv[value+1])
            except:
                print("Syntax Exception: Invalid score value")
    return path, min_loop_len, score_gc, score_au, score_gu

def read_fasta(path):
    sequences = []
    headers = []
    seq = ""
    with open(path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if len(headers) != 0:
                    sequences.append(seq)
                    seq = ""
                headers.append(line)
            else:
                seq = seq + line
        sequences.append(seq)
    return headers, sequences

def fill_matrix(sequence, min_loop, score_gc, score_au, score_gu):
    mat = []
    for i in range(0, len(sequence)):
        mat.append([0]*len(sequence))
    for i in range(1, len(sequence)):
        for j in range(0, len(sequence)-i):
            po_1 = mat[i+1][j]
            po_2 = mat[i][j-1]
            if (sequence[i] == 'A' and sequence[j] == 'U') or (sequence[i] == 'U' and sequence[j] == 'A'):
                po_3 = mat[i+1][j-1] + score_au
            elif (sequence[i] == 'G' and sequence[j] == 'U') or (sequence[i] == 'U' and sequence[j] == 'G'):
                po_3 = mat[i+1][j-1] + score_gu
            elif (sequence[i] == 'G' and sequence[j] == 'C') or (sequence[i] == 'C' and sequence[j] == 'G'):
                po_3 = mat[i+1][j-1] + score_gc
            else:
                po_3 = mat[i+1][j-1] + 0
            po_4 = -1
            for k in range(i, j):
                if mat[i][k] + mat[k+1][j] > po_4:
                    po_4 = mat[i][k] + mat[k+1][j]
            mat[i][j]=max(po_1, po_2, po_3, po_4)
    print(mat)

if __name__ == '__main__':
    path, min_loop_len, score_gc, score_au, score_gu = read_data_in()
    headers, sequences = read_fasta(path)
    fill_matrix(sequences[0], min_loop_len, score_gc, score_au, score_gu)