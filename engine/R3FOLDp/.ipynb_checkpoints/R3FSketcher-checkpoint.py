
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from scipy.optimize import minimize, approx_fprime
from scipy.sparse.csgraph import floyd_warshall
import numpy as np
import networkx as nx

try:
    from R3FUtils import SeqAndDBNToGraph, DBNToPairs, BreaksLine, ModifiedLine
    from R3FUtils import SeqToIdlist, BaseToResidue, AtomAtomDistance, ACGU
    from R3FUtils import BPS, BPBPS, Transform, Kabsch, BACKBONE_ATOMS, RMSD
    from R3FUtils import CoordsToMatrix, DISTS
except:
    from .R3FUtils import SeqAndDBNToGraph, DBNToPairs, BreaksLine, ModifiedLine
    from .R3FUtils import SeqToIdlist, BaseToResidue, AtomAtomDistance, ACGU
    from .R3FUtils import BPS, BPBPS, Transform, Kabsch, BACKBONE_ATOMS, RMSD
    from .R3FUtils import CoordsToMatrix, DISTS


def WingsToCoords(left, right, basesize = 3):
    """left and right wings (two strings) -> stem coords"""

    L = len(left)
    assert L == len(right), "Unequal wings: {} {}".format(left, right)

    ind = L - 1
    bp  = left[ind] + right[ind]
    coords = BPS[bp] if bp in BPS else BPS["GC"] # allow arbitrary base pairs
    nxt    = coords

    while ind > 0:

        ind -= 1
        bp   = left[ind] + right[ind]
        #### GCAU -> APPROPRIATE STEP
        prv_pos  = Transform(BPBPS['GCAU']['PREV'],
                             Kabsch(nxt, BPBPS['GCAU']['NEXT']))
        # GC to allow arbitrary base pairs
        bp_template = BPS[bp] if bp in BPS else BPS["GC"]
        prv = Transform(bp_template,
                        Kabsch(prv_pos, bp_template))

        coords = np.vstack((prv[:basesize],coords,prv[basesize:]))
        nxt    = prv

    return coords
    

def SeqDbnToStems(seq, dbn, basesize, flanked = False):
    """seq+dbn -> paired(set) + stems(list: {'IDLIST':[i,], 'COORDS': Lx3}"""
    
    paired  = set()
    stems   = [] 
    curstem = []

    N = len(seq)
    
    pairs  = DBNToPairs(dbn)
    
    for bp in pairs:

        paired.add(bp[0])
        paired.add(bp[1])
        
        if not curstem:
            curstem.append(bp)
        # if two consecutive base pairs with no chain breaks
        elif bp[0] == curstem[-1][0] + 1 and bp[1] + 1 == curstem[-1][1] and\
             seq[curstem[-1][0]].isupper() and seq[bp[1]].isupper():
            curstem.append(bp)
        else:
            stems.append(curstem)
            curstem = [bp,]

    if curstem:
        stems.append(curstem)

    ##############EXTEND WITH FLANKS#####################
    if flanked:
        for i in range(len(stems)):
            first, last = stems[i][0], stems[i][-1]
            if first[0]-1 >= 0 and first[1]+1 < N:
                stems[i] = [(first[0]-1,first[1]+1),] + stems[i]
            if (last[1]-1) - (last[0]+1) > 2:
                stems[i] = stems[i] + [(last[0]+1,last[1]-1),]
    #####################################################

    for i in range(len(stems)):
        left  = ''.join([seq[bp[0]] for bp in stems[i]]).upper()
        right = ''.join([seq[bp[1]] for bp in stems[i]]).upper()
        #::3 - for C1' only
        coords = WingsToCoords(left, right)[::3] if basesize==1 else WingsToCoords(left, right)
        stems[i] = {'IDLIST': sorted([nt for bp in stems[i] for nt in bp]),
                    'COORDS': coords,
                    'TYPE':   'STEM'}
        

    return stems, paired        


def SeqDbnToFragments(seq, dbn, basesize):
    """split seq+dbn into bases and stems"""
    
    fragments = []

    stems, paired = SeqDbnToStems(seq, dbn, basesize)

    for stem in stems:
        fragments.append(stem)

    for i in range(len(seq)):
        if i not in paired:
            if basesize == 1:
                fragments.append({'IDLIST':[i,],'COORDS': ACGU[seq[i].upper()]['COORDS3'][::3], 'TYPE': 'BASE'})
            else:
                fragments.append({'IDLIST':[i,],'COORDS': ACGU[seq[i].upper()]['COORDS3'], 'TYPE': 'BASE'})

    return fragments


def SeqDbnToApproxFragments(seq, dbn, basesize):
    """split seq+dbn into steps and flanks"""
    
    fragments = []
    N = len(seq)

    stems, paired = SeqDbnToStems(seq, dbn, basesize, flanked = True)

    for stem in stems:
        fragments.append(stem)

    #TODO: GCAU -> sequence-specific steps
    if seq[-1].isupper():
        fragments.append({'IDLIST':[0,N-1],
                          'COORDS': np.vstack((BPBPS['GCAU']['NEXT'][:basesize],
                                               BPBPS['GCAU']['PREV'][:basesize])),
                          'TYPE': 'STEP'})
    for i in range(1,N):
        if seq[i-1].isupper():
            fragments.append({'IDLIST':[i-1, i],
                          'COORDS': np.vstack((BPBPS['GCAU']['PREV'][:basesize],
                                               BPBPS['GCAU']['NEXT'][:basesize])),
                          'TYPE': 'STEP'})

    return fragments  

def GetExactDistances(fragments, restraints, N, BS, MAXW, randomstart, RAND, DIM = 3):
    """fragments   - [list of residues/stems {"IDLIST":..., "COORDS":..., "TYPE": "BASE"/"STEM"}]
       restraints  - [list of restraints {"IDLIST":..., "COORDS":..., "STRUCT": entry}]
       N           - len(seq)
       BS          - basesize
       MAXW        - weight for exact distances
       randomstart - start with random coords (instead of preparing starting coords explicitly)
       RAND        - shifts for starting coords are picked randomly from (-RAND,RAND) range

    """

    M  = np.zeros((BS*N, BS*N)) # Matrix
    W  = np.zeros((BS*N, BS*N)) # Weights  

    if randomstart:
        startcoords = None
    else:
        startcoords = np.zeros((BS*N, DIM))

    #APPLY base+stem fragments to get 0-MAXW weighted 3Nx3N matrix
    for frag in fragments:

        if not randomstart:
            startcoords[[BS*k+b for k in frag['IDLIST']
                                for b in range(BS)]] = frag['COORDS'] - frag['COORDS'][0] + [np.random.randint(-RAND,RAND),
                                                                                             np.random.randint(-RAND,RAND),
                                                                                             np.random.randint(-RAND,RAND)]
        for i in range(len(frag['IDLIST'])):
            k1 = frag['IDLIST'][i]
            for j in range(i,len(frag['IDLIST'])):
                k2 = frag['IDLIST'][j]
                for b1 in range(BS):
                    for b2 in range(BS):
                        D = AtomAtomDistance(frag['COORDS'][BS*i + b1],
                                             frag['COORDS'][BS*j + b2])
                        M[BS*k1 + b1, BS*k2 + b2] = D
                        M[BS*k2 + b2, BS*k1 + b1] = D
                        W[BS*k1 + b1, BS*k2 + b2] = MAXW
                        W[BS*k2 + b2, BS*k1 + b1] = MAXW

    ##### APPLYING RESTRAINTS
    if restraints:
        counts = np.zeros((BS*N, BS*N))
        for frag in restraints:

            if not randomstart:
                startcoords[[BS*k+b for k in frag['IDLIST']
                                    for b in range(BS)]] = frag['COORDS'] - frag['COORDS'][0] + [np.random.randint(-RAND,RAND),
                                                                                                 np.random.randint(-RAND,RAND),
                                                                                                 np.random.randint(-RAND,RAND)]
            for i in range(len(frag['IDLIST'])):
                k1 = frag['IDLIST'][i]
                for j in range(i,len(frag['IDLIST'])):
                    k2 = frag['IDLIST'][j]

                    for b1 in range(BS):
                        for b2 in range(BS):
                            D = AtomAtomDistance(frag['COORDS'][BS*i + b1],
                                                 frag['COORDS'][BS*j + b2])
                            prevcount = counts[BS*k1 + b1, BS*k2 + b2]
                            counts[BS*k1 + b1, BS*k2 + b2] += 1
                            # if this restraint is the first (or unique)
                            if prevcount == 0:
                                M[BS*k1 + b1, BS*k2 + b2] = D
                                M[BS*k2 + b2, BS*k1 + b1] = D
                                W[BS*k1 + b1, BS*k2 + b2] = MAXW
                                W[BS*k2 + b2, BS*k1 + b1] = MAXW
                            else:
                                #if overlapping restraints - do averaging
                                newD = (prevcount*M[BS*k1 + b1, BS*k2 + b2] + D)/(prevcount + 1)
                                M[BS*k1 + b1, BS*k2 + b2] = newD
                                M[BS*k2 + b2, BS*k1 + b1] = newD
    return M, W, startcoords


def GetApproxDistances(M, W, fragments, N, BS, maxw):

    #APPLY step+flank fragments
    for frag in fragments:
        for i in range(len(frag['IDLIST'])):
            k1 = frag['IDLIST'][i]
            for j in range(i,len(frag['IDLIST'])):
                k2 = frag['IDLIST'][j]
                
                if not W[BS*k1, BS*k2] > 0:
                    
                    for b1 in range(BS):
                        for b2 in range(BS):
                            D = AtomAtomDistance(frag['COORDS'][BS*i + b1],
                                                 frag['COORDS'][BS*j + b2])
                            M[BS*k1 + b1, BS*k2 + b2] = D
                            M[BS*k2 + b2, BS*k1 + b1] = D
                            W[BS*k1 + b1, BS*k2 + b2] = maxw
                            W[BS*k2 + b2, BS*k1 + b1] = maxw
    return M, W


def IsMirrored(initcoords, fragments, restraints, basesize):
    
    mirrored = set()
    totalset = set()
    #checking stems & restraints for mirroring
    for frag in fragments +\
                ([] if not restraints
                    else restraints):
        idlist = [basesize*k+b for k in frag['IDLIST'] for b in range(basesize)]
        mirrored_rmsd = RMSD(initcoords[idlist],frag['COORDS']*[-1, 1, 1],superimpose = True)
        rmsd          = RMSD(initcoords[idlist],frag['COORDS'],superimpose = True)
        fragset = set(frag['IDLIST'])
        totalset |= fragset
        if mirrored_rmsd < rmsd:
            mirrored |= fragset
            #print('y',fragset)
            #initcoords[idlist] *= [-1, 1, 1]
        #else:
        #   print('n',fragset)

    #print(len(mirrored)/len(totalset))
    #print(mirrored)
    return len(mirrored) > len(totalset) / 2


def SketchModel(seq, dbn, restraints = None, basesize = 1,
                MAXW = 10, STEP = 3.0, maxw = 1, wpower = 2,
                randomstart = True, RAND = 10, skew = 0.1,
                maxiter = 1, converg = 0.9, dpower = 1,
                mirror = True):
    
    fragments = SeqDbnToFragments(seq, dbn, basesize)

    M, W, startcoords = GetExactDistances(fragments, restraints,
                                          len(seq), basesize,
                                          MAXW, randomstart, RAND)

    ########STEPS AND FLANKS####
    approxfragments = SeqDbnToApproxFragments(seq, dbn, basesize)

    M, W = GetApproxDistances(M, W, approxfragments,
                              len(seq), basesize,
                              maxw)
    ############################

    # DERIVE the weighted matrix of shortest paths
    S = floyd_warshall(M, directed = False)
    M = M*(W == MAXW) + S*(W != MAXW)
    W += 10**-6*(W == 0)
    ##############################################

    return M, W



    
