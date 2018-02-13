# NW.py
# HW1, Computational Genomics, Fall 2018
# andrewid:

# WARNING: Do not change the file name; Autograder expects it.

import sys,math

def ReadFASTA(filename):
    fp=open(filename, 'r')
    Sequences={}
    tmpname=""
    tmpseq=""
    for line in fp:
        if line[0]==">":
            if len(tmpseq)!=0:
                Sequences[tmpname]=tmpseq
            tmpname=line.strip().split()[0][1:]
            tmpseq=""
        else:
            tmpseq+=line.strip()
    Sequences[tmpname]=tmpseq
    fp.close()
    return Sequences

# You may define any helper functions for Needleman-Wunsch algorithm here
def build_matrix(seq1,seq2):
    score = 0
    print_seq1,print_seq2 = '',''
    m = [[0 for row in range(len(seq1)+1)]for col in range(len(seq2)+1)]
    m[0][0] = 0
    i= 0
    while i < len(seq1)+1:
        m[0][i] = -1 * i
        i = i + 1
    j = 1
    while j < len(seq2)+1:
        m[j][0] = -1 * j
        j = j + 1
    #print m

    for index2 in range(1,len(seq2)+1):
        for index1 in range(1,len(seq1)+1):
            if seq1[index1-1] == seq2[index2-1]:
                count_score_nogap = m[index2-1][index1-1]

            elif seq1[index1-1] != seq2[index2-1]:
                count_score_nogap = -1 + m[index2-1][index1-1]
            count_score_gap1 = -1 + m[index2][index1-1]
            count_score_gap2 = -1 + m[index2-1][index1]
            choosemax = max(count_score_nogap,count_score_gap1,count_score_gap2)
            if choosemax == count_score_nogap:
                m[index2][index1] = count_score_nogap
                #print 'loop1'
            elif choosemax == count_score_gap1:
                m[index2][index1] = count_score_gap1
                #print 'loop2'
            elif choosemax == count_score_gap2:
                m[index2][index1] = count_score_gap2
                #print 'loop3'
    return m

def trace_matrix(seq1,seq2,m):
    score = int(m[len(m) - 1][len(m[0]) - 1])
    #print 'score is',score
    x = len(m) - 1
    #print x
    y = len(m[0]) - 1
    #print y
    new_seq1 = ''
    new_seq2 = ''
    while x > 0 or y > 0:

        score_gapx = m[x][y-1]
        score_gapy = m[x-1][y]
        score_nogap = m[x-1][y-1]

        max_score = max(score_nogap,score_gapx,score_gapy)
        #print 'max_score is',max_score
        if max_score == score_nogap:
            new_seq2 = seq2[x - 1] + new_seq2
            new_seq1 = seq1[y - 1] + new_seq1
        #    print 'new_seq is',new_seq1,new_seq2
        #    print x,y
            #if x!= 1 and y != 1:
            y = y - 1
            x = x - 1
        elif max_score == score_gapx:

            new_seq2 = '-' + new_seq2
            new_seq1 = seq1[y - 1] + new_seq1
            #print 'new_seq is',new_seq1,new_seq2
            #print x,y
            #if y != 1:
            y = y - 1
        elif max_score == score_gapy:

            new_seq2 = seq2[x - 1] + new_seq2
            new_seq1 =  '-' + new_seq1
            #print 'new_seq is',new_seq1,new_seq2
            #print x,y
            #if x != 1:
            x = x - 1
    return new_seq1,new_seq2,score
# Do not change this function signature

def needleman_wunsch(seq1, seq2):
    """Find the global alignment for seq1 and seq2
    Returns: a string of three lines like so:
    '<alignment score>\n<alignment in seq1>\n<alignment in seq2>'
    """
    m = build_matrix(seq1,seq2)
    new_seq1,new_seq2,score = trace_matrix(seq1,seq2,m)
    return str(score) + '\n' + new_seq1 + '\n' + new_seq2


if __name__=="__main__":
    Sequences=ReadFASTA(sys.argv[1])
    assert len(Sequences.keys())==2, "fasta file contains more than 2 sequences."
    seq1=Sequences[list(Sequences.keys())[0]]
    seq2=Sequences[list(Sequences.keys())[1]]
    print(needleman_wunsch(seq1, seq2))
