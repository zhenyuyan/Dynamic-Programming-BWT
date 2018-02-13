
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

#print bwt_encode('banana')
#print rotate('abcde#')
print bwt_decode('bnn{aa}a')
