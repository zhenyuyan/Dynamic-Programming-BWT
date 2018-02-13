import re
'''
inputValue=input("please input a int data :")
if type(inputValue)!=type(1):
    raise ValueError
else:
    print inputValue
print max(1,2,3,4)
seq1 = 'ATGC'
seq2 = 'ATTGA'
print '\n'
'''
'''
def rle(s):
    """Run Length Encoder
    Args: s, string to be encoded
    Returns: RLE(s)
    """
    index = 0
    newString = ''
    while index < len(s) - 1:
        if s[index] != s[index + 1] :
            newString = newString + s[index]
            index = index + 1
        elif s[index] == s[index + 1]:
            newString = newString + s[index] + s[index]
            count = 1
            while index < len(s) - 1 and s[index] == s[index + 1]:
                count = count + 1
                index = index + 1
            index = index + 1
            newString = newString + str(count)

    if s == None:
        raise NotImplementedError
    return newString

print rle('WOOOOOHOOOOHOOOO!')
'''
'''
s = 'WOOOOOHOOOOHOOOO!'
s = s.strip(s[0])
print s,s[0],s[1]

def bwt_decode(bwt):
    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """
    if bwt == None:
        raise NotImplementedError
    s1 = bwt.strip('{')
    s2 = s1.strip('}')

    numbers = re.sub("\D", "", s2)
    index = 0
    decoded_string = ''
    num_index = 0
    if len(numbers) == 0:
        return s2
    else:
        while index < len(s2):
            if num_index < len(numbers) and s2[index] != numbers[num_index]:
                decoded_string = decoded_string + s2[index]
            #    print decoded_string
                index = index + 1
                continue
            elif num_index < len(numbers) and s2[index] == numbers[num_index]:
                num = int(numbers[num_index])
                i = 2
                while i < num:
                    decoded_string = decoded_string + s2[index-1]
                    i = i + 1
                if num_index != len(numbers) - 1:
                    num_index = num_index + 1
                index = index + 1
                #print decoded_string
                continue
            index = index + 1
        return decoded_string


s = '{' + 'scott2ytartanscott2ytartanscott2ytartanscott2ytartan' + '}'
print bwt_decode(s)

def count_score(seq1,seq2,score):
    index = 0
    tail_score = abs(len(seq1)-len(seq2))
    while index < min(len(seq1),len(seq2)):
        if seq1[index] != seq2[index]:
            score = score - 1
        index = index + 1
    score = score - tail_score
    return score

def needleman_wunsch(seq1, seq2):
    """Find the global alignment for seq1 and seq2
    Returns: a string of three lines like so:
    '<alignment score>\n<alignment in seq1>\n<alignment in seq2>'
    """
    if seq1 == None or seq2 == None:
        raise NotImplementedError
    print_seq1 = ''
    print_seq2 = ''
    index1 = 0
    index2 = 0
    score = 0
    while index1 < len(seq1):
        while index2 < len(seq2):
            count_score_nogap = 0
            count_score_gap1 = 0
            count_score_gap2 = 0

            if seq1[index1] == seq2[index2]:
                count_score_nogap = 0 + count_score(seq1[index1+1:],seq2[index2+1:],score)
            elif seq1[index1] != seq2[index2]:
                count_score_nogap = -1 + count_score(seq1[index1+1:],seq2[index2+1:],score)

            count_score_gap1 = -1 + count_score(seq1[index1:],seq2[index2+1:],score)
            count_score_gap2 = -1 + count_score(seq1[index1+1:],seq2[index2:],score)

            choosemax = max(count_score_nogap,count_score_gap1,count_score_gap2)

            if choosemax == count_score_nogap and seq1[index1] == seq2[index2]:
                print_seq1 = print_seq1 + seq1[index1]
                print_seq2 = print_seq2 + seq2[index2]
                index1 = index1 + 1
                index2 = index2 + 1

            elif choosemax == count_score_nogap and seq1[index1] != seq2[index2]:
                print_seq1 = print_seq1 + seq1[index1]
                print_seq2 = print_seq2 + seq2[index2]
                score = score - 1
                index1 = index1 + 1
                index2 = index2 + 1

            elif choosemax == count_score_gap1:
                print_seq1 = print_seq1 + '-'
                print_seq2 = print_seq2 + seq2[index2]
                score = score -1
                index2 = index2 + 1

            elif choosemax == count_score_gap2:
                print_seq1 = print_seq1 + seq1[index1]
                print_seq2 = print_seq2 + '-'
                score = score - 1
                index1 = index1 + 1
    newScore = str(score)
    return newScore,print_seq1,print_seq2
'''

seq1 = "TAGCTAATA"
seq2 = "TAGATA"

seq1 = 'ATGC'
seq2 = 'ATTGA'

def build_matrix(seq1, seq2):
    """Find the global alignment for seq1 and seq2
    Returns: a string of three lines like so:
    '<alignment score>\n<alignment in seq1>\n<alignment in seq2>'
    """
    index1 = 1
    index2 = 1
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
    #        print "index1 index2 ",index1,index2
            if seq1[index1-1] == seq2[index2-1]:
                count_score_nogap = m[index2-1][index1-1]
            #    print "no gap equal"
            elif seq1[index1-1] != seq2[index2-1]:
            #    print "no gap not equal"
            #    print 'index not equal',m[index2][index1],index2,index1
                count_score_nogap = -1 + m[index2-1][index1-1]
            count_score_gap1 = -1 + m[index2][index1-1]
#            print 'here',m[index2][index1-1],index2,index1 - 1

            count_score_gap2 = -1 + m[index2-1][index1]
#            print 'ffff',m[index2-1][index1],index2-1,index1

#            print 'count_score_nogap',count_score_nogap
#            print 'count_score_gap1',count_score_gap1
#            print 'count_score_gap2',count_score_gap2
            choosemax = max(count_score_nogap,count_score_gap1,count_score_gap2)
#            print 'max is',choosemax
            if choosemax == count_score_nogap:
                m[index2][index1] = count_score_nogap
#                print 'loop1'
            elif choosemax == count_score_gap1:
                m[index2][index1] = count_score_gap1
#                print 'loop2'
            elif choosemax == count_score_gap2:
                m[index2][index1] = count_score_gap2
#                print 'loop3'
            #m[index1 + 1][index2 + 1] = max(-1 + m[index1][index2 + 1], -1 + m[index1][index2 + 1],-1 + m[index1][index2],m[index1][index2])
    return m

def trace_matrix(seq1,seq2,m):
    score = int(m[len(m) - 1][len(m[0]) - 1])
    print 'score is',score
    x = len(m) - 1
    print x
    y = len(m[0]) - 1
    print y
    new_seq1 = ''
    new_seq2 = ''
    while x > 0 or y > 0:

        score_gapx = m[x][y-1]
        score_gapy = m[x-1][y]
        score_nogap = m[x-1][y-1]

        max_score = max(score_nogap,score_gapx,score_gapy)
        print 'max_score is',max_score
        if max_score == score_nogap:
            new_seq2 = seq2[x - 1] + new_seq2
            new_seq1 = seq1[y - 1] + new_seq1
            print 'new_seq is',new_seq1,new_seq2
            print x,y
            #if x!= 1 and y != 1:
            y = y - 1
            x = x - 1
        elif max_score == score_gapx:

            new_seq2 = '-' + new_seq2
            new_seq1 = seq1[y - 1] + new_seq1
            print 'new_seq is',new_seq1,new_seq2
            print x,y
            #if y != 1:
            y = y - 1
        elif max_score == score_gapy:

            new_seq2 = seq2[x - 1] + new_seq2
            new_seq1 =  '-' + new_seq1
            print 'new_seq is',new_seq1,new_seq2
            print x,y
            #if x != 1:
            x = x - 1


    return new_seq1,new_seq2

m = build_matrix(seq1,seq2)
print m
print (trace_matrix(seq1,seq2,m))
'''
       T    A   G   C   T   A   A   T   A
[ [-0, -1, -2, -3, -4, -5, -6, -7, -8, -9],
T [-1, -0, -1, -2, -3, -4, -5, -6, -7, -8],
A [-2, -1, -0, -1, -2, -3, -4, -5, -6, -7],
G [-3, -2, -1, -0, -1, -2, -3, -4, -5, -6],
A [-4, -3, -2, -1, -1, -2, -2, -3, -4, -5],
T [-5, -4, -3, -2, -2, -1, -2, -3, -3, -4],
A [-6, -5, -4, -3, -3, -2, -1, -2, -3, -3]]
#print (needleman_wunsch(seq1,seq2))
'''
