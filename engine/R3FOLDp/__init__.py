
import numpy as np
import os

try:
    from .R3FUtils  import CoordsToMatrix, PairsToDBN, BASE3, Base
    from .R3FSketcher import SketchModel
    from .R3FParser import FileToEntry, EntryToStructure, SplitStructures
except:
    from R3FUtils  import CoordsToMatrix, PairsToDBN, BASE3, Base
    from R3FSketcher import SketchModel
    from R3FParser import FileToEntry, EntryToStructure, SplitStructures


def ParseOneAtom(file):

    res = []

    entry = FileToEntry(file, fileformat = 'pdb')

    #IGNORE EMPTY FILES
    if len(entry['IDLIST']) < 7:
        return res
    ###################
    
    structs = EntryToStructure(entry)
    for struct in SplitStructures(structs):

        coords = np.array([struct['IDDICT'][strid][BASE3[Base(strid)][0]]
                           for strid in struct['IDLIST']])
        MY = CoordsToMatrix(coords)

        #REMOVING LONE BASE PAIRS
        #pairs = struct['BPS']
        #nolone = []
        #for i in range(len(pairs)):
        #    if i > 0 and (pairs[i-1][0]+1 == pairs[i][0] and\
        #                  pairs[i-1][1]-1 == pairs[i][1]) or\
        #       i < len(pairs) - 1 and (pairs[i][0]+1 == pairs[i+1][0] and\
        #                               pairs[i][1]-1 == pairs[i+1][1]):
        #        nolone.append(pairs[i])
        #########################

        seq = struct['SEQ']
        dbn = PairsToDBN(struct['BPS'],len(seq))
        MX, MW = SketchModel(seq,dbn)
        
        res.append([seq,dbn,struct['BPS'],MX,MW,MY,os.path.basename(file)[:-4]])

    return res
        
