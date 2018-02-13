# BWT.py
# HW1, Computational Genomics, Fall 2018
# andrewid:zhenyuy1

# WARNING: Do not change the file name, or the function signatures below.
# Autograder expects these names exactly.
import re,sys
def rle(s):
    """Run Length Encoder
    Args: s, string to be encoded
    Returns: RLE(s)
    """
    if s == None:
        raise NotImplementedError
    s = s.strip('{').strip('}')
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
    if s[len(s)- 1] != s[len(s) - 2]:
        newString = newString + s[len(s)-1]
    return newString

def rotate(s):
    i = 0
    l = [s]
    while i < len(s) - 1:
       s = s[len(s) - 1] + s[0:len(s)-1]
       l.append(s)
       i = i + 1
    return l

def bwt_encode(s):
    """Burrows-Wheeler Transform
    Args: s, string, which must not contain '{' or '}'
    Returns: BWT(s), which contains '{' and '}'
    """
    if s == None:
        raise NotImplementedError
    newString = '{' + s + '}'

    l = rotate(newString)
    l.sort()
    encoded_word = ''
    for word in l:
        encoded_word = encoded_word + word[len(word) - 1]
    return encoded_word


def bwt_decode(bwt):
    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """
    if bwt == None:
        raise NotImplementedError
    #bwt = bwt.strip('{').strip('}').strip('#')
    #bwt = bwt + '#'
    #print bwt
    l0 = list(bwt)
    l = []
    for word in l0:
       l.append(word)

    for i in range(len(bwt) - 1):
         l.sort()
         j = 0
         while j < len(l):
             l[j] = l0[j] + l[j]
             j = j + 1
    choose_word = ''
    for word in l:
        if word[0] == '{':
            choose_word = word
    choose_word = choose_word.strip('{').strip('}')
    return choose_word

# test_strings is not graded, just provided for your own convenience
def test_string(s):
    compressed = rle(s)
    bwt = bwt_encode(s)
    compressed_bwt = rle(bwt)
    reconstructed = bwt_decode(bwt)
    template = "{:25} ({:3d}) {}"
    print(template.format("original", len(s), s))
    print(template.format("bwt_enc(orig)", len(bwt), bwt))
    print(template.format("bwt_dec(bwt_enc(orig))", len(reconstructed), reconstructed))
    print(template.format("rle(orig)", len(compressed), compressed))
    print(template.format("rle(bwt_enc(orig))", len(compressed_bwt), compressed_bwt))
    print()
    print()

if __name__ == "__main__":
    # Add more of your own strings to explore for question (i)
    test_strings = [sys.argv[1],
                    sys.argv[2]]
    for s in test_strings:
        test_string(s)
