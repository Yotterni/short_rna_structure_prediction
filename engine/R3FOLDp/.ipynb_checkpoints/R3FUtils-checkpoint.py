
import numpy as np

# Gapped values
GAPS = {'-', '.', '~'}

BACKBONE_ATOMS = set("P OP1 OP2 C1' C2' O2' C3' O3' C4' O4' C5' O5'".split())


BASE_SRC = """
    
A      C1'    C2    C6    A
G      C1'    C2    C6    G

1MA    C1'   C2    C6    A
2MG    C1'   C2    C6    G
6MZ    C1'   C2    C6    A
7MG    C1'   C2    C6    G
G7M    C1'   C2    C6    G
A2M    C1'   C2    C6    A
MA6    C1'   C2    C6    A
OMG    C1'   C2    C6    G
YYG    C1'   C2    C6    G
SAM    C1'   C2    C6    A

C      C1'   C2    C4    C
U      C1'   C2    C4    U

4AC    C1'   C2    C4    C
4SU    C1'   C2    C4    U
5MC    C1'   C2    C4    C
5MU    C1'   C2    C4    U
H2U    C1'   C2    C4    U
LV2    C1'   C2    C4    C
OMC    C1'   C2    C4    C
OMU    C1'   C2    C4    U
SSU    C1'   C2    C4    U
UR3    C1'   C2    C4    U

PSU    C1'   C4    C2    U
B8N    C1'   C4    C2    U
3TD    C1'   C4    C2    U
UY1    C1'   C4    C2    U

DA     C1'   C2    C6    A
DG     C1'   C2    C6    G
DC     C1'   C2    C4    C
DT     C1'   C2    C4    U
DU     C1'   C2    C4    U
    
ADP    C1'   C2    C6    A
GDP    C1'   C2    C6    G
CDP    C1'   C2    C4    C
UDP    C1'   C2    C4    U

ATP    C1'   C2    C6    A
GTP    C1'   C2    C6    G
CTP    C1'   C2    C4    C
UTP    C1'   C2    C4    U

12A	C1'	C2	C6	A
16B	C1'	C2	C4	C
1MA	C1'	C2	C6	A
1MG	C1'	C2	C6	G
1RN	C1'	C2	C4	U
1RP	C1'	CAQ	NAL	U
1SC	C1'	C2	C4	C
1W5	C1'	C2	C4	C
23G	C1'	C2	C6	G
2IA	C1'	C2	C6	A
2MA	C1'	C2	C6	A
2MG	C1'	C2	C6	G
2MU	C1'	C2	C4	U
2PR	C1'	C2	C6	G
2SG	C1'	C2	C6	G
3AU	C1'	C2	C4	U
3DA	C1'	C2	C6	A
3TD	C1'	C4	C2	U
4AC	C1'	C2	C4	C
4DU	C1'	C2	C6	A
4OC	C1'	C2	C4	C
4SU	C1'	C2	C4	U
56B	C1'	C2	C6	G
5BU	C1'	C2	C4	U
5CF	C1'	C2	C4	C
5FU	C1'	C2	C4	U
5GP	C1'	C2	C6	G
5HM	C1'	C2	C4	C
5IC	C1'	C2	C4	C
5MC	C1'	C2	C4	C
5MU	C1'	C2	C4	U
5UD	C1'	C2	C4	U
6FC	C1'	C2	C4	C
6IA	C1'	C2	C6	A
6MD	C1'	C2	C6	A
6MZ	C1'	C2	C6	A
70U	C1'	C2	C4	U
73W	C1'	C2	C4	C
75B	C1'	C2	C4	U
7AT	C1'	C2	C6	A
7MG	C1'	C2	C6	G
7S3	C1'	C2	C6	G
7SN	C1'	C2	C6	G
8AZ	C1'	C2	C6	G
8B4	C1'	C2	C6	A
8OS	C6	C2	C1	G
8RJ	C1G	C2	C4	U
9QV	C1'	C2	C4	U
A	C1'	C2	C6	A
A23	C1'	C2	C6	A
A2M	C1'	C2	C6	A
A5M	C1'	C2	C4	C
A7C	C1'	C2	C6	A
ADP	C1'	C2	C6	A
ADS	C1'	C2	C6	A
AET	C1'	C2	C6	A
AF2	C1'	C2	C6	A
AG9	C1'	C2	C4	C
AMP	C1'	C2	C6	A
AMZ	C1	C7A	N1	C
APC	C1'	C2	C6	A
AT9	C1'	C2	C6	G
ATP	C1'	C2	C6	A
AVC	C1'	C2	C6	A
B8H	C1'	C4	C2	U
B8N	C1'	C4	C2	U
B8T	C1'	C2	C4	C
BGM	C1'	C2	C6	G
BRU	C1'	C2	C4	U
C	C1'	C2	C4	C
C4J	C1'	C4	C2	U
C5P	C1'	C2	C4	C
CBR	C1'	C2	C4	C
CBV	C1'	C2	C4	C
CCC	C1'	C2	C4	C
CFL	C1'	C2	C4	C
CFZ	C1'	C2	C4	C
CH1	C1'	C2	C4	C
CM0	C1'	C2	C4	U
CSL	C1'	C2	C4	C
CTP	C1'	C2	C4	C
D4M	C1'	C2	C4	U
DA	C1'	C2	C6	A
DC	C1'	C2	C4	C
DG	C1'	C2	C6	G
DGP	C1'	C2	C6	G
DI	C1'	C2	C6	G
DOC	C1'	C2	C4	C
DSH	C1'	C2	C6	A
DT	C1'	C2	C4	U
DU	C1'	C2	C4	U
F86	C1	C12	C11	A
FHU	C1'	C4	C2	U
FMU	C1'	C2	C4	U
G	C1'	C2	C6	G
G2L	C1'	C2	C6	G
G46	C1'	C2	C6	G
G47	C1'	C2	C6	G
G4P	C1'	C2	C6	G
G5J	C1'	C2	C6	G
G7M	C1'	C2	C6	G
GAO	C1'	C2	C6	G
GDO	C1'	C2	C6	G
GDP	C1'	C2	C6	G
GE6	C03	C10	N03	C
GF2	C1'	C2	C6	G
GH3	C1'	C2	C6	G
GMP	C1'	C2	C6	G
GMX	C1'	C2	C6	G
GNG	C1'	C2	C6	G
GRB	C1'	C2	C6	G
GTP	C1'	C2	C6	G
H2U	C1'	C2	C4	U
HCU	C6	C2	C1	G
HYJ	C1'	C2	C6	G
I	C1'	C2	C6	G
I2T	C1'	C4	C2	U
IC5	C1'	C2	C6	G
IKS	C1'	C2	C4	C
ILK	C1'	C2	C4	U
IU	C1'	C2	C4	U
JMC	C1'	C2	C4	C
JMH	C1'	C2	C4	C
K1F	C1'	C22	N05	C
L2B	C1'	C2	C4	U
L3X	C1'	C2	C6	A
LCA	C1'	C2	C6	A
LCC	C1'	C2	C4	C
LCG	C1'	C2	C6	G
LHH	C1'	C2	C4	C
LKC	C1'	C2	C4	C
LV2	C1'	C2	C4	C
LXR	C5	C6	C9	A
M1Y	C1'	C4	C2	U
M2G	C1'	C2	C6	G
M3X	C1'	C2	C4	C
M7A	C1'	C2	C6	A
M7G	C1'	C2	C6	G
M7M	CBQ	CBM	CBF	G
MA6	C1'	C2	C6	A
MGT	C1'	C2	C6	G
MIA	C1'	C2	C6	A
MMX	C1'	C2	C4	C
MUM	C1'	C2	C4	U
NF2	C1'	C2	C4	U
NMN	C1R	C2	C4	C
O2C	C1'	C2	C4	C
OMC	C1'	C2	C4	C
OMG	C1'	C2	C6	G
OMU	C1'	C2	C4	U
P5P	C1'	C2	C6	A
PGP	C1'	C2	C6	G
PPU	C1'	C2	C6	A
PSU	C1'	C4	C2	U
PYO	C1'	C2	C4	U
PYY	C1'	C2	C4	U
PZG	C6	C2	C1	G
QUO	C1'	C2	C6	G
RIA	C1A	C2	C6	A
RPC	C1'	C2	C4	C
RSP	C1'	C2	C4	C
RSQ	C1'	C2	C4	C
RUS	C1'	C2	C4	U
RY	C1'	C2	C4	C
S4C	C1'	C2	C4	C
SAH	C1'	C2	C6	A
SAM	C1'	C2	C6	A
SUR	C1'	C2	C4	U
T6A	C1'	C2	C6	A
TLN	C1'	C2	C4	U
TPP	C7'	C4'	C2'	C
TTP	C1'	C2	C4	U
U	C1'	C2	C4	U
U23	C1'	C2	C4	U
U33	C1'	C2	C4	U
U5M	C1'	C2	C4	U
U5P	C1'	C2	C4	U
U8U	C1'	C2	C4	U
UBD	C1'	C2	C4	U
UD5	C1'	C2	C4	U
UFP	C1'	C2	C4	U
UFT	C1'	C2	C4	U
UMS	C1'	C2	C4	U
UR3	C1'	C2	C4	U
US5	C1'	C2	C4	U
UTP	C1'	C2	C4	U
UVP	C1'	C2	C4	U
UY1	C1'	C4	C2	U
UY4	C1'	C2	C6	A
W0F	C06	C15	C17	G
XMP	C1'	C2	C6	G
XUG	C1'	C2	C6	G
YG	C1'	C2	C6	G
YYG	C1'	C2	C6	G
ZJS	C1'	C2	C6	A
""".strip().split('\n')


#BASE3 = DICT {'A': ['N9', 'C2', 'C6'], ...}
BASE3    = {x.split()[0]:x.split()[1:4] for x in BASE_SRC if x.strip()}
#MODIFIED = DICT {'A': 'A', 'SAM': 'A', ...}
MODIFIED = {x.split()[0]:x.split()[4]   for x in BASE_SRC if x.strip()}


ACGU_SRC = """
P   A           -8.195  20.363  15.990
OP1 A           -8.086  21.020  17.307
OP2 A           -6.993  19.992  15.206
O5' A           -9.054  19.055  16.176
C5' A           -9.979  18.896  17.282
C4' A          -10.409  17.427  17.268
O4' A          -10.928  17.123  15.935
C3' A           -9.377  16.365  17.455
O3' A           -8.970  16.320  18.826
C2' A          -10.160  15.187  17.014
O2' A          -11.201  14.848  17.890
C1' A          -10.714  15.745  15.706
N9  A           -9.793  15.604  14.592
C8  A           -8.893  16.492  14.063
N7  A           -8.234  15.988  13.040
C5  A           -8.731  14.699  12.893
C6  A           -8.455  13.644  12.003
N6  A           -7.560  13.687  11.016
N1  A           -9.160  12.503  12.172
C2  A          -10.071  12.402  13.145
N3  A          -10.413  13.326  14.038 
C4  A           -9.696  14.446  13.844 
P   C            0.253   8.558   8.361 
OP1 C            0.690   9.131   9.639 
OP2 C           -0.311   7.190   8.365 
O5' C           -0.815   9.491   7.657 
C5' C           -0.548  10.875   7.448 
C4' C           -1.659  11.427   6.579
O4' C           -1.672  10.807   5.284
C3' C           -3.064  11.195   7.081
O3' C           -3.349  12.180   8.078
C2' C           -3.903  11.330   5.814
O2' C           -4.083  12.661   5.395 
C1' C           -3.018  10.630   4.834
N1  C           -3.240   9.193   4.588
C2  C           -4.148   8.833   3.589
O2  C           -4.848   9.711   3.073 
N3  C           -4.257   7.521   3.258 
C4  C           -3.506   6.598   3.856 
N4  C           -3.648   5.324   3.466 
C5  C           -2.620   6.931   4.912
C6  C           -2.524   8.236   5.216
P   G          -10.278   1.916  -4.832 
OP1 G          -10.830   1.788  -6.160
OP2 G           -8.934   1.413  -4.561 
O5' G          -10.247   3.456  -4.374
C5' G          -11.403   4.255  -4.589
C4' G          -11.130   5.623  -4.008
O4' G          -10.993   5.504  -2.589 
C3' G           -9.852   6.307  -4.437 
O3' G          -10.034   6.874  -5.726 
C2' G           -9.624   7.299  -3.314 
O2' G          -10.483   8.406  -3.519
C1' G          -10.052   6.460  -2.133 
N9  G           -8.921   5.762  -1.502 
C8  G           -8.511   4.461  -1.694 
N7  G           -7.500   4.121  -0.978 
C5  G           -7.192   5.264  -0.240 
C6  G           -6.196   5.514   0.744
O6  G           -5.359   4.739   1.202 
N1  G           -6.215   6.843   1.199 
C2  G           -7.089   7.797   0.754
N2  G           -6.903   9.045   1.235
N3  G           -8.047   7.576  -0.142
C4  G           -8.045   6.280  -0.574 
P   U           -8.653   4.653   6.379 
OP1 U           -9.451   3.621   5.688
OP2 U           -8.237   5.906   5.679
O5' U           -9.449   5.159   7.631
C5' U          -10.151   4.162   8.430
C4' U          -10.608   4.897   9.700
O4' U           -9.479   5.661  10.206
C3' U          -11.675   5.928   9.574
O3' U          -12.951   5.271   9.475
C2' U          -11.502   6.681  10.831
O2' U          -12.022   5.976  11.935
C1' U           -9.989   6.769  10.935
N1  U           -9.387   7.976  10.358
C2  U           -9.468   9.107  11.151 
O2  U          -10.004   9.104  12.246 
N3  U           -8.898  10.222  10.597
C4  U           -8.271  10.333   9.374 
O4  U           -7.813  11.426   9.037 
C5  U           -8.230   9.121   8.619 
C6  U           -8.776   8.006   9.126
""".strip().split('\n')

ACGU = {}
for x in ACGU_SRC:
    atomname, base, X, Y, Z = x.strip().split()
    coords = [float(X), float(Y), float(Z)]
    if base not in ACGU:
        ACGU[base] = {'ATOMS': [],'COORDS': {}}
    ACGU[base]['COORDS'][atomname] = coords
    ACGU[base]['ATOMS'].append(atomname)
for base in ACGU:
    ACGU[base]['COORDS3'] = np.array([ACGU[base]['COORDS'][atom]
                                      for atom in BASE3[base]])
    ACGU[base]['COORDS']  = np.array([ACGU[base]['COORDS'][atom]
                                      for atom in ACGU[base]['ATOMS']])


BPS = {'GC': np.array([[-10.052,  6.460,  -2.133],
                       [-7.089,   7.797,   0.754],
                       [-6.196,   5.514,   0.744],
                       [-3.018,  10.630,   4.834],
                       [-4.148,   8.833,   3.589],
                       [-3.506,   6.598,   3.856]]),
       'CG': np.array([[-3.018,  10.630,   4.834],
                       [-4.148,   8.833,   3.589],
                       [-3.506,   6.598,   3.856],
                       [-10.052,  6.460,  -2.133],
                       [-7.089,   7.797,   0.754],
                       [-6.196,   5.514,   0.744]]),
       'AU': np.array([[-10.714, 15.745,  15.706],
                       [-10.071, 12.402,  13.145],
                       [-8.455,  13.644,  12.003],
                       [-9.989,   6.769,  10.935],
                       [-9.468,   9.107,  11.151],
                       [-8.271,  10.333,   9.374]]),
       'UA': np.array([[-9.989,   6.769,  10.935],
                       [-9.468,   9.107,  11.151],
                       [-8.271,  10.333,   9.374],
                       [-10.714, 15.745,  15.706],
                       [-10.071, 12.402,  13.145],
                       [-8.455,  13.644,  12.003],]),
       'GU': np.array([[16.116,  13.540, -15.957],
                       [14.246,  17.265, -17.257],
                       [12.083,  16.098, -17.070],
                       [10.304,  21.875, -18.509],
                       [10.164,  19.472, -18.181],
                       [8.098,   18.226, -18.503]]),
       'UG': np.array([[10.304,  21.875, -18.509],
                       [10.164,  19.472, -18.181],
                       [8.098,   18.226, -18.503],
                       [16.116,  13.540, -15.957],
                       [14.246,  17.265, -17.257],
                       [12.083,  16.098, -17.070],]),}


BPBPS = {'GCAU':{'PREV': np.array([[-4.133,   4.840,  -2.468],
                                   [-1.033,   1.842,  -2.225],
                                   [-0.415,   2.727,  -0.056],
                                   [3.000,  -2.890,  -2.858],
                                   [1.935,   -0.976,  -1.859],
                                   [2.816,    0.071,   0.034]]),
                 'NEXT': np.array([[-4.563,  -0.587,  -2.328],
                                   [-0.982,  -1.924,  -0.250],
                                   [-1.306,  -0.229,   1.317],
                                   [3.907,  -5.136,   1.956],
                                   [2.264,   -3.313,   1.664],
                                   [2.338,   -1.340,   3.062]]),}
         }

def AtomAtomDistance(coords1, coords2):

    return pow(sum((c1 - c2)**2
                   for c1, c2 in zip(coords1, coords2)), 0.5)


def Base(idstr : '1.B.G.25.'):
    """
    Parse the base from a dssr-like id string
    """
    return idstr.split('.')[2]


def Kabsch(arr1, arr2):
    """Return translation & rotation matrices for
       superimposing arr2 to arr1"""

    mean1 = arr1.mean(axis = 0)
    mean2 = arr2.mean(axis = 0)

    X = arr1 - mean1
    Y = arr2 - mean2

    S, _, D = np.linalg.svd(np.dot(np.transpose(Y), X))
    rot = np.transpose(np.dot(np.transpose(D), np.transpose(S)))
    if np.linalg.det(rot) < 0:
        D[2] = -D[2]
        rot = np.transpose(np.dot(np.transpose(D), np.transpose(S)))

    return mean1, mean2, rot
    

def Transform(arr, operators):
    """Apply translation and rotation to arr"""
    mean1, mean2, rot = operators
    return np.dot(arr - mean2, rot) + mean1
    

def RMSD(arr1, arr2, superimpose = False):
    """RMSD between two matrices, if superimpose = True,
       the second matrix will be superimposed to the first one
       before calculating RMSD"""
    if superimpose:
        arr2 = Transform(arr2, Kabsch(arr1, arr2))
    dX = arr1 - arr2    
    return np.sqrt(np.mean(np.sum(np.multiply(dX, dX),axis=1)))


def ModelModelTransform(model1, model2):
    """Superimpose model2 to model1"""
    return {'SEQ': model2['SEQ'], 'DBN': model2['DBN'],
            'COORDS': Transform(model2['COORDS'],
                                Kabsch(model1['COORDS'], model2['COORDS']))}
    

def ModelModelRMSD(model1, model2, superimpose = True):
    """3-atom RMSD between two models (by default, superimpose first)"""
    return RMSD(model1['COORDS'],
                model2['COORDS'],
                superimpose = superimpose)


def StructStructRMSD(struct1, struct2, atom = "C3'", superimpose = True):
    """1-atom (by default C3') RMSD between two structures
       (by default, superimpose first)"""
    arr1 = []
    arr2 = []

    for idstr1, idstr2 in zip(struct1['IDLIST'], struct2['IDLIST']):
        arr1.append(struct1['IDDICT'][idstr1][atom])
        arr2.append(struct2['IDDICT'][idstr2][atom])

    return RMSD(np.array(arr1),
                np.array(arr2),
                superimpose = superimpose)


def BasePairReferenceDistance(bp, ref):
    """return two distances:
       RMSD from one two-residue superposition &
       RMSD from two one-residue superpositions"""
    refbp = BPS[ref]

    rmsd1 = RMSD(refbp, bp, True)
    rmsd2 = (RMSD(refbp[3:], Transform(bp[3:], Kabsch(refbp[:3],bp[:3]))) +\
             RMSD(refbp[:3], Transform(bp[:3], Kabsch(refbp[3:],bp[3:])))) / 2

    return rmsd1, rmsd2


def BreaksLine(seq):
    """list of binary 3'-ending/non-3'-ending residue values"""
    return [int(x.islower()) for x in seq]


def CleanIdList(entry):
    """removes non-residue idstrs from idlist
    and returns the cleaned idlist
    (residues with missing B1/B2/B3 atoms
     are treated as non-residues)"""
    res = []

    for idstr in entry['IDLIST']:

        base = Base(idstr)
        if base in BASE3 and all(atom in entry['IDDICT'][idstr]
                                 for atom in BASE3[base]):
            res.append(idstr)
    
    return res


def IntToChain(n):
    """Ordered chain identifiers: A,B,...,z,AA,AB,...,zz"""
    chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz'

    if n < 0:
        raise Exception("ERROR: negative input for IntToChain func")

    if n < len(chars):
        return ' '+chars[n]

    if n > len(chars)**2 + len(chars) - 1:
        raise Exception("ERROR: IntToChain cannot handle value {} > {}"\
                        .format(n,len(chars)**2 + len(chars) - 1))

    return chars[(n // len(chars)) - 1] + chars[n % len(chars)]


def ModifiedLine(idlist):
    """list of binary modified/standard residue values"""
    line = []
    for idstr in idlist:
        base = Base(idstr)
        line.append(1 - int(base in MODIFIED and MODIFIED[base] == base))
    return line


def DBNToPairs(dbn):
    """Convert the dbn string into a sorted list of base pairs"""
    pairs = set()

    # keys == closing brackets, values == matching opening brackets
    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D',
               'e':'E','f':'F','g':'G','h':'H','i':'I','j':'J','k':'K','l':'L',
               'm':'M','n':'N','o':'O','p':'P','q':'Q','r':'R','s':'S','t':'T',
               'u':'U','v':'V','w':'W','x':'X','y':'Y','z':'Z',
               'б':'Б','г':'Г','д':'Д','ё':'Ё','ж':'Ж','й':'Й','л':'Л','п':'П',
               'ф':'Ф','ц':'Ц','ч':'Ч','ш':'Ш','щ':'Щ','ь':'Ь','ы':'Ы','ъ':'Ъ',
               'э':'Э','ю':'Ю','я':'Я',}
    # 30+19 bp stacks for 30+19 allowed pseudoknot levels
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[],
             'F':[],'G':[],'H':[],'I':[],'J':[],'K':[],'L':[],'M':[],'N':[],
             'O':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'U':[],'V':[],'W':[],
             'X':[],'Y':[],'Z':[],
             'Б':[],'Г':[],'Д':[],'Ё':[],'Ж':[],'Й':[],'Л':[],'П':[],'Ф':[],
             'Ц':[],'Ч':[],'Ш':[],'Щ':[],'Ь':[],'Ы':[],'Ъ':[],'Э':[],'Ю':[],
             'Я':[],}
              
    for i,v in enumerate(dbn):
        # if we observe an opening bracket
        # then add its index into the matching stack
        if v in stack: 
            stack[v].append(i)
        # else if we observe the closing bracket
        # take the opening index from the matching stack
        # and add the base pair to the pairs set
        elif v in closing:
            # this is to handle closing brackets with no
            # opening partner - they will be ignored
            if stack[closing[v]]:
                pairs.add((stack[closing[v]].pop(), i))

    return sorted(pairs)
      

def PairsToDBN(newpairs, length = 0, returnlevels = False, levellimit = -1):
    """Convert a list of base pairs into a dbn string of the given length"""

    # Initialize the dbn string
    dbn = ['.']*length

    # Define "brackets" for 30 pseudoknot levels (and 19 more encoded with cyrillic letters)
    # Higher levels will be simply ignored
    levels = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff','Gg',
              'Hh','Ii','Jj','Kk','Ll','Mm','Nn','Oo','Pp','Qq','Rr',
              'Ss','Tt','Uu','Vv','Ww','Xx','Yy','Zz',
              'Бб','Гг','Дд','Ёё','Жж','Йй','Лл','Пп',
              'Фф','Цц','Чч','Шш','Щщ','Ьь','Ыы','Ъъ','Ээ','Юю','Яя']

    # groups of non-conflicting base pairs
    groups = [set(),]

    # normalize the pairs (i.e. ensure v < w)
    pairs = set((min(v, w), max(v, w)) for v, w in newpairs)
    
    for pair in sorted(pairs):

        level = 0

        # find the minimum level where the pair is not in conflict
        # with any base pair of that level
        while any(v[0]<=pair[0]<=v[1]<=pair[1] or
                  pair[0]<=v[0]<=pair[1]<=v[1] for v in groups[level]):
            level += 1
            if level == len(groups):
                groups.append(set())
            if level == len(levels):
                levels.append('..')

        # add the pair to the determined level
        groups[level].add(pair)

    # kind of a bubble sort of the base pairs among the levels
    # to maximize the number of base pairs of the lowest levels
    # e.g. to turn (..[[[...)...]]] into [..(((...]...)))
    for times in range(len(groups)-1):
        for i in range(len(groups)-1):

            rest = {v for v in groups[i+1] if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                       w[0]<=v[0]<=w[1]<=v[1]
                                                       for w in groups[i])}
            clean = groups[i+1] - rest

            while rest:

                confjprev = set()
                confiprev = set()

                confj = rest.pop()
                rest.add(confj)
                confj = {confj,}
                confi = {v for v in groups[i] if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                     w[0]<=v[0]<=w[1]<=v[1]
                                                     for w in confj)}

                while confjprev != confj or confiprev != confi:

                    confjprev = confj
                    confiprev = confi

                    confj = {v for v in rest if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                    w[0]<=v[0]<=w[1]<=v[1]
                                                    for w in confi)}
                    confi = {v for v in groups[i] if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                     w[0]<=v[0]<=w[1]<=v[1]
                                                     for w in confj)}

                if len(confi) < len(confj):

                    groups[i]   = confj | (groups[i] - confi)
                    groups[i+1] = confi | (groups[i+1] - confj)

                rest = rest - confj

            if clean:

                groups[i] |= clean
                groups[i+1] -= clean

    if returnlevels:
        levels = {}
        for lev, group in enumerate(groups):
            for bp in group:
                levels[bp] = lev + 1
        return levels

    # remove all levels higher than levellimit (if specified)
    if levellimit >= 0:
        groups = groups[:levellimit]

    # add all the pairs to the dbn string
    # according to their levels  
    for i, group in enumerate(groups):
        for pair in group:
            dbn[pair[0]] = levels[i][0]
            dbn[pair[1]] = levels[i][1]
            
    return ''.join(dbn)


def AdjToCC(adj):
    """adjacency list to connected components (bfs)"""
    ccs   = []
    seen  = set()
    queue = []

    while len(seen) < len(adj):
        cur   = set()
        queue.append([k for k in adj.keys() if k not in seen][0])
        while queue:
            v = queue[0]
            queue = queue[1:]
            if v not in seen:
                seen.add(v)
                cur.add(v)
                for w in adj[v]:
                    if w not in seen:
                        queue.append(w)
        ccs.append(cur)
    return ccs

def SeqAndDBNToGraph(seq, dbn):
    """ seq & dbn -> {'ADJ': adjacency list (dictionary) {i: {j,...},},
                      'CC': connected components [{i,j,...},{...}]}"""
    adj = {i:{} for i in range(len(seq))}
    ccs = []

    N = len(seq)

    #covalent bonds
    for i in range(len(seq)):
        if not seq[i - 1].islower():
            v = i
            # w = previous or last in case of 0th
            w = (len(seq) - 1) if v == 0 else (v - 1)
            adj[v][w] = 1.0
            adj[w][v] = 1.0

    for v, w in DBNToPairs(dbn):
        adj[v][w] = 2.0
        adj[w][v] = 2.0

        '''if v+1 < w and seq[v].isupper():
            adj[v+1][w] = 1.5
            adj[w][v+1] = 1.5
            if v+2 < w and seq[v+1].isupper():
                adj[v+2][w] = 1.5
                adj[w][v+2] = 1.5
        if w+1 < N and seq[w].isupper():
            adj[v][w+1] = 1.5
            adj[w+1][v] = 1.5
            if w+2 < N and seq[w+1].isupper():
                adj[v][w+2] = 1.5
                adj[w+2][v] = 1.5'''

    ccs = AdjToCC(adj)

    return {"ADJ": adj, "CC": ccs}


def SeqToIdlist(seq):
    """"Gg" -> ["1.A.G.1.", "1.A.G.2."]"""
    cnt_asym  = 0
    cnt_seq   = 1
    cur_chain = IntToChain(cnt_asym).strip()

    idlist = []

    for base in seq:

        idlist.append('.'.join(['1',          #Model
                                cur_chain,    #Chain
                                base.upper(), #Base
                                str(cnt_seq), #Residue id
                                '']))         #Ins.code  
        cnt_seq += 1

        if base.islower():
            cnt_asym += 1
            cur_chain = IntToChain(cnt_asym).strip()
            cnt_seq = 1

    return idlist


def BaseToResidue(idstr, coords, raw = False, backbone = None):
    """idstr + 3-atom coords -> full-atom residue coords"""
    base = Base(idstr)

    residue = {}

    #in case we want original 3-atom coords in the coordinate file
    if raw:
        residue['ATOMS'] = [x for x in BASE3[base]]
        for i in range(len(residue['ATOMS'])):
            residue[residue['ATOMS'][i]] = list(coords[i])
        return residue
    #####

    residue['ATOMS'] = [x for x in ACGU[base]['ATOMS']]

    rb_atoms_d = {}

    if backbone:
        ratomdict, ridstr = backbone
        r3_atoms = [x for x in BASE3[Base(ridstr)]]
        rb_atoms = [x for x in ratomdict if x in BACKBONE_ATOMS]
        rb_atoms_d = {x:i for i,x in enumerate(rb_atoms)}

        r3_coords = np.array([ratomdict[x] for x in r3_atoms])
        rb_coords = np.array([ratomdict[x] for x in rb_atoms])

        backbone_coords = Transform(rb_coords,
                                    Kabsch(coords, r3_coords))

        transformed_r3_coords = Transform(r3_coords,
                                          Kabsch(coords, r3_coords))

        # shift to fix N1/N9 positions exactly
        res3_coords = Transform(ACGU[base]['COORDS3'],
                                Kabsch(coords, ACGU[base]['COORDS3']))
        shift = res3_coords[0] - transformed_r3_coords[0]
        backbone_coords += shift


    res_coords = Transform(ACGU[base]['COORDS'],
                           Kabsch(coords, ACGU[base]['COORDS3']))

    for i in range(len(residue['ATOMS'])):

        atomname = residue['ATOMS'][i]

        if atomname in rb_atoms_d:
            residue[atomname] = list(backbone_coords[rb_atoms_d[atomname]])
        else:
            residue[atomname] = list(res_coords[i])

    return residue


def CoordsToMatrix(coords):
    """3*N array -> 3*Nx3*N distance matrix"""
    return np.sum((coords[:, np.newaxis, :] - coords[np.newaxis, :, :]) ** 2,
                  axis = -1)**0.5


def IsSeqConnected(seq, dbn):
    """seq & dbn -> bool connected/not connected"""
    return len(SeqAndDBNToGraph(seq, dbn)['CC']) == 1


DISTS_SRC = """1.0 6.07989211161019 4604
1.5 9.27386645016948 2748
2.0 10.610670204806596 5495
2.5 12.592628206994677 3302
3.0 14.458430290762593 6042
3.5 15.792280697623719 3585
4.0 17.46659855445336 6507
4.5 18.657613171410894 4136
5.0 19.765087876492398 6786
5.5 21.355204523515464 4740
6.0 21.669066803843887 7048
6.5 24.00715752135093 5181
7.0 23.306664822975897 7338
7.5 26.18411116918818 5778
8.0 24.97116886653913 7608
8.5 28.10320300701169 6366
9.0 26.797192268076156 7864
9.5 29.89319056894681 6980
10.0 28.748965842674814 8156
10.5 31.67409155443813 7593
11.0 30.818853005096347 8443
11.5 33.48006178898882 8059
12.0 32.973020266422715 8740
12.5 35.217799810719676 8542
13.0 35.08603868487282 9028
13.5 36.87932577069283 9036
14.0 36.94482614511493 9407
14.5 38.518288291936166 9583
15.0 38.67721048394572 9806
15.5 40.1267837192041 10099
16.0 40.297660711592094 10218
16.5 41.588667697011275 10579
17.0 41.81090560066232 10630
17.5 43.021401253443 11001
18.0 43.25454955834223 10914
18.5 44.38415981424564 11284
19.0 44.63219027872781 11215
19.5 45.704463468926505 11611
20.0 45.96141684511722 11588
20.5 46.999146699796995 11983
21.0 47.28821899260973 11901
21.5 48.318874997877344 12282
22.0 48.68095547723265 12168
22.5 49.623295806886155 12593
23.0 50.01740802485842 12516
23.5 50.89080317203058 12976
24.0 51.35533968570491 12951
24.5 52.257742104360105 13367
25.0 52.6346858650368 13380
25.5 53.627050097544426 13805
26.0 53.969647173281714 13804
26.5 54.87691064438718 14286
27.0 55.18444392879587 14359
27.5 56.10213811388409 14789
28.0 56.36808361832199 14980
28.5 57.357975176948514 15331
29.0 57.57807333188129 15595
29.5 58.593053337313776 15932
30.0 58.70020148277662 16222
30.5 59.76000165317643 16567
31.0 59.81391858012056 16861
31.5 60.78099579585081 17198
32.0 60.930871471405425 17412
32.5 61.79633997843806 17736
33.0 62.01442073457653 17954
33.5 62.89305052494134 18334
34.0 63.035493628677 18431
34.5 63.91805263583536 18921
35.0 64.04282938687366 18998
35.5 64.9579515849403 19499
36.0 65.14074078892521 19582
36.5 65.9862499372153 20125
37.0 66.18765233099771 20082
37.5 66.92903169298863 20625
38.0 67.16725916985956 20481
38.5 67.91420633801509 21053
39.0 68.26149861245266 20894
39.5 68.8457390938185 21488
40.0 69.15133433815954 21261
40.5 69.77379997561344 21879
41.0 70.11916235360313 21615
41.5 70.67238077695791 22243
42.0 71.0761238621668 21943
42.5 71.54550084785653 22586
43.0 71.97519889023522 22265
43.5 72.34588240647251 22859
44.0 72.8970502370152 22646
44.5 73.16064757750593 23226
45.0 73.72499473001758 23093
45.5 73.91331432066983 23682
46.0 74.56394318917923 23596
46.5 74.63796463098808 24118
47.0 75.30741857991023 24112
47.5 75.33415388772535 24576
48.0 76.09715806101073 24527
48.5 76.08919642015525 25039
49.0 76.85531569272473 24982
49.5 76.84920034059469 25512
50.0 77.61325528049696 25427
50.5 77.57817793296051 26037
51.0 78.3808658379635 25923
51.5 78.3838778211808 26516
52.0 79.10592031501409 26502
52.5 79.11970775651383 27054
53.0 79.88927579458979 27127
53.5 79.87649763905699 27806
54.0 80.56562207360241 27777
54.5 80.66955191459877 28435
55.0 81.31763103953966 28469
55.5 81.43103637198246 29006
56.0 82.06294257141612 28997
56.5 82.08481459984377 29581
57.0 82.78900549537909 29495
57.5 82.68442895210039 30199
58.0 83.51193418122564 30040
58.5 83.29580215309211 30794
59.0 84.19879627209049 30624
59.5 83.96442652729301 31368
60.0 84.8366766284375 31069
60.5 84.6745392027558 31720
61.0 85.47721242834878 31382
61.5 85.35724630405493 31962
62.0 86.08810103353638 31658
62.5 85.94906365723794 32097
63.0 86.67146894144646 31790
63.5 86.43966811556659 32210
64.0 87.22629802852066 31912
64.5 86.93819137922789 32300
65.0 87.74000823100205 31922
65.5 87.35013502687481 32318
66.0 88.1834683469335 31887
66.5 87.77557238858624 32229
67.0 88.68047060498891 31773
67.5 88.21360213350019 32040
68.0 89.19251445284162 31564
68.5 88.7434403882485 31830
69.0 89.67972938636554 31400
69.5 89.2817581436087 31642
70.0 90.22690249890411 31210
70.5 89.81913657377318 31426
71.0 90.79471929865785 30947
71.5 90.33616571338626 31138
72.0 91.32066348025178 30684
72.5 90.91348789083843 30889
73.0 91.74081342884145 30515
73.5 91.42217797501596 30569
74.0 92.10190719108626 30313
74.5 91.91330892523035 30212
75.0 92.51816197731556 30049
75.5 92.36464433379403 29912
76.0 92.86389856951112 29785
76.5 92.80584455191924 29676
77.0 93.26883059473015 29544
77.5 93.1564607246859 29410
78.0 93.60597689450687 29321
78.5 93.41432290076575 29222
79.0 93.83126052489605 29114
79.5 93.64111829362332 29145
80.0 94.00135353645061 28944
80.5 93.84553080905647 28933
81.0 94.13083039835966 28687
81.5 94.01748187764116 28733
82.0 94.32786433316029 28374
82.5 94.12819951315997 28464
83.0 94.53017568129134 28063
83.5 94.23892406886148 28231
84.0 94.76267098841224 27789
84.5 94.45209736603371 27960
85.0 94.9648556548128 27582
85.5 94.7284357927081 27622
86.0 95.22004341852673 27368
86.5 95.06531876547288 27222
87.0 95.52244499288157 27091
87.5 95.40984130236754 26866
88.0 95.80835904382998 26781
88.5 95.80461397177146 26463
89.0 96.14550992233059 26422
89.5 96.22420648018593 26063
90.0 96.51741423316724 25977
90.5 96.60541465744163 25656
91.0 96.87357975621987 25524
91.5 97.02269501459462 25264
92.0 97.35457448474341 25114
92.5 97.48342275631674 24866
93.0 97.76073815356408 24687
93.5 97.94445251894032 24336
94.0 98.17154268780467 24198
94.5 98.51090242599902 23708
95.0 98.73174756607132 23621
95.5 99.18251619776748 23006
96.0 99.21731944935829 22960
96.5 99.74488516840033 22263
97.0 99.68405025357387 22275
97.5 100.28054272013856 21603
98.0 100.17544247301954 21671
98.5 100.97161213092699 20967
99.0 100.72164067358467 21053
99.5 101.51833402920188 20298
100.0 101.1163412011992 20436
100.5 101.9582935235176 19682
101.0 101.42156963961594 19839
101.5 102.4377193340349 19017
102.0 101.77618806918476 19279
102.5 102.9637311375314 18449
103.0 102.0040456623634 18732
103.5 103.43869467342478 17889
104.0 102.43135180001985 18138
104.5 103.92928888901295 17413
105.0 102.77259223603154 17611
105.5 104.40543923803605 16996
106.0 103.10734355827442 17132
106.5 104.85819180298249 16565
107.0 103.46985083729835 16623
107.5 105.37242111849537 16185
108.0 103.89324551018245 16164
108.5 105.85551908986386 15839
109.0 104.218553128515 15766
109.5 106.41095012079106 15477
110.0 104.65756219101831 15360
110.5 106.96722709711447 15104
111.0 105.15732065491437 14995
111.5 107.54413285801385 14753
112.0 105.62363553785313 14631
112.5 108.07719390350302 14389
113.0 106.18070335991564 14280
113.5 108.58902459711253 13986
114.0 106.70881844897612 13876
114.5 109.10372886799757 13584
115.0 107.36356854485061 13458
115.5 109.51596621458086 13202
116.0 107.99908397887592 13006
116.5 110.11591025637154 12784
117.0 108.58109709719875 12587
117.5 110.69091277919225 12373
118.0 109.19353513847999 12176
118.5 111.3559660451048 11956
119.0 109.88974109449201 11744
119.5 112.08703542476323 11538
120.0 110.69910579127357 11304
120.5 112.83613953164736 11091
121.0 111.40976265538534 10891
121.5 113.37301131185804 10677
122.0 111.96015445653997 10521
122.5 113.8727196282946 10292
123.0 112.46791901471472 10178
123.5 114.32810619350103 9901
124.0 112.96922160747661 9794
124.5 114.67578692510565 9535
125.0 113.48618696539175 9418
125.5 115.13461539139709 9198
126.0 114.0296694138868 9048
126.5 115.56012617890445 8877
127.0 114.66014555534355 8724
127.5 116.07085086335125 8572
128.0 115.43238795225056 8401
128.5 116.7191952099222 8264
129.0 116.25759292314453 8117
129.5 117.21025753892131 7986
130.0 116.90351774041447 7845
130.5 117.85304530091005 7698
131.0 117.48895571653392 7588
131.5 118.14873184598027 7415
132.0 117.90685980419008 7319
132.5 118.51995297178152 7164
133.0 118.45206810685269 7052
133.5 119.01696013883918 6932
134.0 119.25939958333865 6788
134.5 119.62971623390999 6698
135.0 120.04377212665666 6537
135.5 120.30175704897707 6468
136.0 120.89641852941277 6286
136.5 120.99632534528301 6222
137.0 121.52906282425583 6072
137.5 121.68036148971272 5945
138.0 122.25027152621944 5815
138.5 122.26910748454397 5622
139.0 123.00424947555356 5543
139.5 122.94828586521396 5303
140.0 123.8255052347322 5284
140.5 123.60129425765501 4998
141.0 124.35826456490233 5025
141.5 123.99274998823184 4706
142.0 124.83542902263973 4773
142.5 124.3717748503431 4432
143.0 125.30469037644143 4516
143.5 124.76395349097214 4168
144.0 125.74186787622038 4247
144.5 124.94187094090056 3904
145.0 125.94445700740765 3997
145.5 125.16315112911917 3645
146.0 126.28186338604351 3736
146.5 125.3814521409118 3388
147.0 126.27687082554455 3522
147.5 125.39871524159662 3172
148.0 126.529467824949 3311
148.5 125.49097355322033 2980
149.0 126.7901508808572 3123
149.5 125.31325050295948 2814
150.0 126.9841163202915 2957
150.5 125.43953291253268 2641
151.0 127.2004438215518 2777
151.5 125.21174809221628 2479
152.0 127.47356328770078 2591
152.5 125.2761392516532 2308
153.0 127.95890408020298 2397
153.5 125.51366269304476 2150
154.0 128.0003923326351 2226
154.5 125.7022022818536 2009
155.0 128.3437965341351 2083
155.5 126.23789051146193 1879
156.0 128.79603853398012 1964
156.5 127.20927992718319 1756
157.0 129.46791747678643 1841
157.5 128.48697510782685 1644
158.0 129.96183816341787 1710
158.5 129.95529272954803 1518
159.0 130.92693502934276 1575
159.5 131.62077116633338 1395
160.0 132.45966733104444 1453
160.5 132.78882990692966 1283
161.0 133.82705478168904 1340
161.5 133.74398338660504 1195
162.0 134.8683199813512 1247
162.5 134.64303697138388 1112
163.0 135.844666421059 1149
163.5 135.37405136923365 1050
164.0 136.8945848148272 1077
164.5 136.47502602764766 1007
165.0 137.11312969487483 1027
165.5 137.5764388694963 973
166.0 137.66416912643118 979
166.5 138.4291944674466 928
167.0 138.1055533362049 936
167.5 140.04928811400995 884
168.0 138.79476842337883 892
168.5 142.04176143930337 836
169.0 139.93908100691044 843
169.5 143.78681622487804 789
170.0 141.0815163477797 791
170.5 145.14910255664287 748
171.0 143.1806203542887 731
171.5 146.58515311779175 701
172.0 145.4839877836262 673
172.5 147.98751647732396 650
173.0 147.07738627010323 609
173.5 148.8006174841669 600
174.0 148.17610270408713 556
174.5 149.60312841078286 552
175.0 148.36087817169687 511
175.5 150.9910257712457 516
176.0 148.4713703950009 466
176.5 152.13016109946838 478
177.0 149.34614134812674 423
177.5 153.58852904183632 438
178.0 149.51922467291106 389
178.5 153.2717006243938 407
179.0 149.5193017521124 360
179.5 153.1994368140404 376
180.0 149.14354703378024 335""".split('\n')

DISTS = {float(x.split()[0]):float(x.split()[1]) for x in DISTS_SRC}


if __name__ == "__main__":

    seq = 'GGGAGcGCGAAGAAcAAaCCcGGg'
    dbn = '(...(())..[...)..(.].)..'
    #seq = 'a'
    #dbn = '.'

    coords = np.array([[0,0],
                       [1,1],
                       [2,2],
                       [2,0]])
    print(CoordsToMatrix(coords))


    #idlist = SeqToIdlist(seq)
    #for x in idlist:
    #    print(x)

    '''g = SeqAndDBNToGraph(seq, dbn)
    print(seq)
    print(dbn)
    print(g['ADJ'])
    print(g['CC'])
    print(IsSeqConnected(seq, dbn))
    print(BASE3)'''

    '''X = np.array([[1,3],[1,5],[3,3]])
    Y = np.array([[4,5],[6,5],[4,3]])
    Z = np.array([[5,4],[7,4],[5,2]])
    print(RMSD(X,Y))
    print(RMSD(X,Y,True))
    print(RMSD(X, Transform(Z, Kabsch(X, Y))))
    print(Transform(Z, Kabsch(X, Y)))
    print(Transform(Y, Kabsch(X, Y)))'''
