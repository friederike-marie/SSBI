"""Dominik Molitor, Friederike Moroff, 05.05.2022, Assignment 01"""
import sys


def read_data_in():
    """Reads in Command-line inputs as options for Nussinov algorithm

    :return str path: Path to file as string
    :return int min_loop_len: lower boundary for loop size
    :return int score_gc: value of GC bonds in Nussinov calculation
    :return int score_au: value of AU bonds in Nussinov calculation
    :return int score_gu: value of non-canonical GU bonds in Nussinov calculation
    """
    # First example FastA file settings as a default
    min_loop_len = 10
    score_gc = 4
    score_au = 1
    score_gu = 0

    for value in range(0, len(sys.argv)):  # For each argument in command line
        if "-i" not in sys.argv:
            print("Error: Path is missing, please enter a correct FastA file")
            quit()
        if sys.argv[value] == "-i":  # input of the file
            if not sys.argv[value + 1].startswith("-"):
                path = sys.argv[value + 1]
            else:  # short check that the file is not empty
                print("Syntax Error: Incorrect syntax at -i")
                quit()
        elif sys.argv[value] == "--min-loop-length":  # if minimal loop length has been entered...
            try:
                min_loop_len = int(sys.argv[value + 1])  # ...read in next argument as valid integer...
            except:  # ...or throw an Exception
                print("Syntax Exception: Invalid minimum loop length\nThe standard length of 3 will be used")
        elif sys.argv[value] == "--score-GC":  # if score GC is entered...
            try:
                score_gc = int(sys.argv[value + 1])   # ...read next argument in as valid integer...
            except:
                print("Syntax Exception: Invalid score value")   # ...or throw an Exception
        elif sys.argv[value] == "--score-AU":
            try:
                score_au = int(sys.argv[value + 1])
            except:
                print("Syntax Exception: Invalid score value")
        elif sys.argv[value] == "--score-GU":  # if score GC is entered...
            try:
                score_gu = int(sys.argv[value + 1])   # ...read next argument in as valid integer...
            except:
                print("Syntax Exception: Invalid score value")   # ...or throw an Exception
    return path, min_loop_len, score_gc, score_au, score_gu


def read_fasta(path):
    """Reads in fasta as headers and sequences

    :param str path: String which acts as the URL

    :return list headers: list of all fasta headers
    :return list sequences: list of all fasta sequences
    """
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
    """Fills a matrix according to a modified Nussinov algorithm

    :param str sequence: RNA sequence on which Nussinov will be performed
    :param int min_loop: limit how many residues have to be at least in a loop ==> distance to create a bond
    :param int score_gc: value of gc bonds, for unmodified Nussinov set to 1
    :param int score_au: value of au bonds, for unmodified Nussinov set to 1
    :param int score_gu: value of gu bonds, for unmodified Nussinov set to 0

    :return list mat: a list of list as a matrix implementation. Every entry is filled or 0
    """
    # Initialisation of the matrix
    mat = []
    for i in range(0, len(sequence)):
        mat.append([0] * len(sequence))

    # Nussinov, for each column run down a diagonal and fill with largest entry
    for column in range(min_loop + 1, len(sequence), 1):
        for row in range(0, len(sequence) - column):
            opt1 = mat[row + 1][column + row]   #option 1: entry below current entry --> no bond
            opt2 = mat[row][column + row - 1]   #option 2: entry left of current entry --> no bond
            # option 3: calculate the score-value + down-left entry in matrix
            if (sequence[row] == "A" and sequence[row + column] == "U") or (sequence[row] == "U" and sequence[
                row + column] == "A"):
                opt3 = mat[row + 1][column + row - 1] + score_au
            elif (sequence[row] == "G" and sequence[row + column] == "C") or (sequence[row] == "C" and sequence[
                row + column] == "G"):
                opt3 = mat[row + 1][column + row - 1] + score_gc
            elif (sequence[row] == "G" and sequence[row + column] == "U") or (sequence[row] == "U" and sequence[
                row + column] == "G"):
                opt3 = mat[row + 1][column + row - 1] + score_gu
            else:
                opt3 = mat[row + 1][column + row - 1] # and if no match just the value
            # option 4: Combine two sub structures with each other with the highest combined value
            opt4 = -1
            for k in range(row + 1, column + row - 1):
                if mat[row][k] + mat[k + 1][row + column] > opt4:
                    opt4 = mat[row][k] + mat[k + 1][row + column]
            # Finally, find maximum
            value = max(opt1, opt2, opt3, opt4)
            mat[row][column + row] = value
    return mat


def create_traceback(matrix, row=0, column=0):
    """Generates a traceback of the Nussinov algorithm

    :param list matrix: list of list filled by fill_matrix
    :param int row: Index in which row the traceback currently is
    :param int column: Index in which column the traceback currently is

    :return list alignment: list of tuples of all bonds
    """

    # Initialisations
    alignment = []
    done = False

    # Traceback
    while matrix[row][column] != 0: # As long as there is atleast one more bond <--> value != 0
        if matrix[row][column] == matrix[row][column - 1]: # If the value left to it is the same move there
            column = column - 1
        elif matrix[row][column] == matrix[row + 1][column]: # If the value below is the same move there
            row = row + 1
        # Check for a bond and move diagonally down-left...
        elif matrix[row][column] == matrix[row + 1][column - 1] + score_au or matrix[row][column] == matrix[row + 1][
            column - 1] + score_gc or matrix[row][column] == matrix[row + 1][column - 1] + score_gu:
            alignment.append((row + 1, column + 1)) # ...add bond to alignment list
            row = row + 1
            column = column - 1
        # If there were two sub-structures continue the traceback of these two substructures
        else:
            for k in range(row + 1, column - 1):
                if matrix[row][column] == matrix[row][k] + matrix[k + 1][column] and not done:
                    alignment += create_traceback(matrix, row=row, column=k)
                    alignment += create_traceback(matrix, row=k + 1, column=column)
                    return alignment
    return alignment


def write_file(alignment, file_name, min_l, gc, au, gu, matrix, seq):
    """Write Output file in a .txt-format

    :param list alignment: list of tuples of the coordinates of bonds
    :param string file_name: Original file name
    :param int min_l: minimum number of residues in loop
    :param int gc: value of gc bonds, for unmodified Nussinov set to 1
    :param int au: value of au bonds, for unmodified Nussinov set to 1
    :param int gu: value of gu bonds, for unmodified Nussinov set to 1
    :param list matrix: List of list which is used as a matrix for Nussinov
    :param str seq: Original RNA sequence
    """
    with open("OutPut.txt", 'a') as file:
        dot_bracket = ""
        file.writelines("Filename: {}\n".format(file_name))
        file.writelines("Min-loop: {}\n".format(min_l))
        file.writelines("GC: {}\n".format(gc))
        file.writelines("AU: {}\n".format(au))
        file.writelines("GU: {}\n".format(gu))
        file.writelines("score: {}\n".format(matrix[0][-1]))
        for i in range(0, len(seq)):
            match = False
            for j in alignment:
                if j[0] == i + 1:
                    dot_bracket = dot_bracket + "("
                    file.writelines(str(j[0]) + "\t" + seq[i] + "\t" + str(j[1]) + "\n")
                    match = True
                elif j[1] == i + 1:
                    dot_bracket = dot_bracket + ")"
                    file.writelines(str(j[1]) + "\t" + seq[i] + "\t" + str(j[0]) + "\n")
                    match = True
            if not match:
                file.writelines(str(i + 1) + "\t" + seq[i] + "\t" + "0" + "\n")
                dot_bracket = dot_bracket + "."
        file.writelines(dot_bracket)

def print_matrix(sequences, matrix):
    """Prints matrix into command line

    :param str sequences: Original RNA sequence
    :param list matrix: List of list --> Nussinov matrix
    """
    first_line = ""
    for i in sequences:
        first_line += "\t" + i
    print(first_line)
    for val, i in enumerate(matrix):
        row = sequences[val]
        for j in i:
            row += "\t" + str(j)
        row += "\n"
        print(row)

#Basic main function
if __name__ == '__main__':
    path, min_loop_len, score_gc, score_au, score_gu = read_data_in()
    headers, sequences = read_fasta(path)
    mst = fill_matrix(sequences[0], min_loop_len, score_gc, score_au, score_gu)
    alignment = create_traceback(mst, column=len(sequences[0]) - 1)
    write_file(alignment, path, min_loop_len, score_gc, score_au, score_gu, mst, sequences[0])
    print_matrix(sequences[0], mst)