from Bio import pairwise2
from Bio import SeqIO
from Bio.Align import substitution_matrices
import numpy as np
import csv
import matplotlib.pyplot as plt
import math

"""
key ideas:
function that returns matrix for each PAMn
function to get global alignment score for each PAMn
matplotlib to plot score
"""

np.set_printoptions(suppress=True)

aaFrequency = {0: 0.087, 1: 0.041, 2: 0.040, 3: 0.047, 4: 0.033, 5: 0.038, 6: 0.050, 7: 0.089, 8: 0.034,
               9: 0.037,
               10: 0.085, 11: 0.081, 12: 0.015, 13: 0.040, 14: 0.051, 15: 0.070, 16: 0.058, 17: 0.010,
               18: 0.030,
               19: 0.065}


def read_file(csvfilename):  # function to read file containing pam1 mutation probability matrix
    rows = []

    with open(csvfilename, encoding='utf-8') as csvfile:
        file_reader = csv.reader(csvfile)
        for row in file_reader:
            rows.append(row)
    return rows


def normalizeCol(M):  # function to normalise column of pam1 matrix to ensure that each row sums to 1
    col_sums = M.sum(axis=0)
    return M / col_sums


def getPAMScoringMatrix(n):  # function to get scoring matrix
    matrix = []
    smallConstant = 0.00000000000001
    data = read_file("pam1.csv")[1:]  # exclude header of file
    for row in data:
        row = row[0].split("\t")
        temp = []
        for i in row:
            temp.append(float(i))
        matrix.append(temp)

    pam1Matrix = normalizeCol(np.array(matrix))  # converting list of list into numpy array type and normalising
    pamNMatrix = normalizeCol(np.linalg.matrix_power(pam1Matrix, n))  # normalising array after n matrix multiplications

    scoringMatrix = []

    for i in range(20):
        newRow = []
        for j in range(20):
            if n != 1:
                ratio = pamNMatrix[i][j] / aaFrequency[i]
                newRow.append(round(10 * math.log10(ratio)))
            else:  # if n is 1, adding a small constant to avoid exception raised by log function
                ratio = pamNMatrix[i][j] / aaFrequency[i] + smallConstant
                newRow.append(round(10 * math.log10(ratio)))
        scoringMatrix.append(newRow)
        # print(list(scoringMatrix[i]))

    return scoringMatrix


def getListOfScores(seq1, seq2, n):
    try:  # if seq parsed in contains more than 1 sequences, use parse() function
        seq1 = SeqIO.read(seq1, "fasta")
    except ValueError:
        seq1 = next(SeqIO.parse(seq1, "fasta"))
    try:
        seq2 = SeqIO.read(seq2, "fasta")
    except ValueError:
        seq2 = next(SeqIO.parse(seq2, "fasta"))

    listOfScores = []
    for i in range(1, n + 1):
        print(f"Calculating Needle score for iteration {i} for {seq1.id} and {seq2.id}")
        scoringMatrix = np.array(getPAMScoringMatrix(i))  # get scoring matrix
        # converting matrix into substitution_matrices array type
        scoringMatrix = substitution_matrices.Array(alphabet="ARNDCQEGHILKMFPSTWYV", dims=2, data=scoringMatrix)
        gapPenalty = np.amin(scoringMatrix)  # getting smallest value in the PAM matrix at that moment to be gap penalty
        gapExtensionPenalty = gapPenalty / 10
        # getting global alignment using substitution matrix
        alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, scoringMatrix, gapPenalty, gapExtensionPenalty)
        score = alignments[0][2]  # indexing score from alignments variable
        listOfScores.append(score)
    return listOfScores, seq1.id, seq2.id


def plotData(seq1, seq2, n, ignoreNegatives=None, plot=None, colour=None):
    """
    :param seq1: Name of FASTA file of first sequence
    :param seq2: Name of FASTA file of second sequence
    :param n: Number of PAM Matrices
    :param ignoreNegatives: Input True if graph plotted ignores negative scores
    :param plot: Input "scatter" or "line" to specify type of graph
    :param colour: Colour of line of graph
    :return: Returns None and plots a graph
    """
    lst, seq1id, seq2id = getListOfScores(seq1, seq2, n)
    axis = [i for i in range(1, n + 1)]
    x = np.array(axis)
    y = np.array(lst)
    if ignoreNegatives:
        x = x[y > 0.0]
        y = y[y > 0.0]
    if plot == "scatter":
        plt.scatter(x, y)
    else:
        plt.setp(plt.plot(x, y), color=colour, linewidth=2.0)  # plotting score against PAM matrix
    plt.title(f"Needle score of {seq1id} and {seq2id} for every PAM matrix")
    plt.xlabel("PAM matrix")
    plt.ylabel("Needleman-Wunsch alignment score")
    plt.show()


# you can try with ignoreNegatives = False and plot = "scatter"
plotData("humanalphaglobin.faa", "bovinealphaglobin.faa", 5, colour='r', ignoreNegatives=False, plot="line")
plotData("humanalphaglobin.faa", "humanmyoglobin.faa", 5, colour='b', ignoreNegatives=True, plot="line")