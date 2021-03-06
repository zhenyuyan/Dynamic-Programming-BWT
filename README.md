# Dynamic-Programming-BWT

1 NW.py:
This file is a program of needleman - wunsch dynamic programming algorithm. You can view the introduction here : https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm
Scores :mismatch: - 1, match 0, gap  - 1.

i. input: a fasta file which contains two strings of DNA. You can input example.fasta here.

ii. output: the result of dynamic programming alignment result(score, seqence1,sequence2).
example input: python NW.py example.fasta
example output:
-8
GGTAGCTTAAAGAAGCAGCATAGGTTTTAGGTGATCCCTCAGCTTAACACAAGGGGAAAATACTTTATAGGCTGGTTTGCAAACTATCATTTGC--TGTTTAGTCAAGGCTGCCAAGAAAACTGTTGGAATTCTTA-GTAATTAGCTCAGCAGCTTGGCTTGAATTCAAAATACCAGCTCTGAAGGGATCCCCATATCAGCTCCAGCGCCAAAAAGCACCTCACTTTAGGGACACACTGAGTATCTTAGAGATTGTTTTCCTCTGTTCTCCCAGGCCTATGACATGGAGCACACTTTCTACAGCAATGGAGAGAAGAAGAAGATTTACATGGAAATTGATCCTGTGACCAGAACTGAAATATTCAGAAGCGGAAATGGCACTGATGAAACATTGGAAGTGCACGACTTTAAAAACGGATACACTGGCATCTACTTCGTGGGTCTTCAAAAATGTTTTATCAAAACTCAGATTAAAGTGATTCCTGAATTTTCTGAACCAGAAGAGGAAATAGATGAGAATGAAGAAATTACCACAACTTTCTTTG
GGTAGCTTAAAGAAGCAGCATAGGTTTTAGGTGATCCCTCAGCTTAACACAAGGGGAAAATACTTTATAGGCTGGTTTGCAAACTATCATTTGCGCTGTTTAGTCAAGGCTGCCAAGAAAACTGTTGGAATTCTTAAGTAATTAGCTCAGCAGCTTGGCTTGAATTCAAAATACCAGCTCTGAAGGGATCC--ATATCAGCTCCAGCGCCAAAAAGCACCTCACTTTAGGGACAC-CTGAGTATCTTAGAGATTGTTTTCCTCTGTTCTCCCAGGCCTATGACATGGAGCACACTTTCTACAGCAATGGAGAGAAGAAGAAGATTTACATGGAAATTGATCCTGTGACCAGAACTGAAATATTCAGAAGCGGAAATGGCACTGATGAAACATTGGAAGTGCACGACTTTAAA--CGGATACACTGGCATCTACTTCGTGGGTCTTCAAAAATGTTTTATCAAAACTCAGATTAAAGTGATTCCTGAATTTTCTGAACCAGAAGAGGAAATAGATGAGAATGAAGAAATTACCACAACTTTCTTTG

2 BWT.py:
Burrows-Wheeler Transform and Run Length Encoding (RLE).
You can consult https://en.wikipedia.org/wiki/Burrows–Wheeler_transform for more details.

i. input: two strings that you want to encode,decode with.

ii. output: The result of rle, BWT-encoding , BWT-decoding.
example input: python BWT.py WOOOOHOOOOHOOOO! scottytartanscottytartan
example output:
original                  ( 16) WOOOOHOOOOHOOOO!
bwt_enc(orig)             ( 18) OOOOOOOOOOOOHHW{}!
bwt_dec(bwt_enc(orig))    ( 16) WOOOOHOOOOHOOOO!
rle(orig)                 ( 13) WOO4HOO4HOO4!
rle(bwt_enc(orig))        ( 11) OO12HH2W{}!


original                  ( 24) scottytartanscottytartan
bwt_enc(orig)             ( 26) ttttssaaccaa{nrryyootttt}n
bwt_dec(bwt_enc(orig))    ( 24) scottytartanscottytartan
rle(orig)                 ( 26) scott2ytartanscott2ytartan
rle(bwt_enc(orig))        ( 31) tt4ss2aa2cc2aa2{nrr2yy2oo2tt4}n

