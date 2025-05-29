
import glob, os
from multiprocessing import Pool

from R3FOLDp import ParseOneAtom


if __name__ == "__main__":

    loc = "../RFdiffusion/compile_all/*/*.pdb"
    files = sorted(glob.glob(loc))

    threads = os.cpu_count()

    lens = {}

    with Pool(threads) as pool:
        for entry in pool.imap(ParseOneAtom, files):
            for data in entry:
                
                seq, dbn, bps, mx, mw, my, name = data
                #RNA sequence, lowercase chars for 3'-ends (if single sequence,
                #                                           only the last char
                #                                           is lowercase)
                # for now we can put all to UPPERCASE
                print(seq)

                #dot-bracket line, just to look at
                print(dbn)

                #sorted list of base pairs (i,j), where i > j (0-based)
                print(bps)

                #input matrix of pairwise distances
                print(mx)

                #weight matrix:
                #10     - confident distances;
                # 1     - intermediate;
                #10**-6 - FW-derived shortest paths
                print(mw)

                #output matrix of real pairwise distances
                print(my)

# lengths distrib:
'''
0 7
7 6
8 7
9 10
10 54
11 112
12 190
13 243
14 208
15 220
16 455
17 142
18 200
19 141
20 240
21 405
22 328
23 254
24 311
25 645
26 592
27 594
28 1019
29 740
30 2517
31 808
32 1058
33 602
34 909
35 661
36 1306
37 901
38 971
39 1009
40 8181
41 848
42 955
43 768
44 1368
45 5279
46 997
47 813
48 1314
49 921
50 262619
51 845
52 1761
53 870
54 782
55 984
56 991
57 813
58 682
59 2074
60 32313
61 704
62 984
63 2188
64 887
65 736
66 789
67 744
68 2005
69 715
70 13405
71 674
72 849
73 1682
74 1136
75 5928
76 2053
77 732
78 676
79 771
80 2133
81 722
82 669
83 741
84 978
85 754
86 611
87 503
88 608
89 822
90 1803
91 569
92 512
93 1420
94 565
95 365
96 1554
97 712
98 677
99 1162
100 15051
116 990
128 990
134 1005
137 885
140 985
166 960
174 1005
175 965
201 660
204 1015
206 945
239 5
240 22360
'''
