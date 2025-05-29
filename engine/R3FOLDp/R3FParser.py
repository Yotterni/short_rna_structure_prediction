
import numpy as np

try:
    from R3FUtils  import BreaksLine, CleanIdList, ModifiedLine, SeqAndDBNToGraph, PairsToDBN
    from R3FUtils  import Base, BASE3, MODIFIED, AtomAtomDistance, BPS, BasePairReferenceDistance
    from R3FUtils  import CoordsToMatrix, GAPS
    from R3FReader import GuessFormat, ParsePDB, ParseCIF
    from R3FWriter import WriteToPDB, WriteToCIF
except:
    from .R3FUtils  import BreaksLine, CleanIdList, ModifiedLine, SeqAndDBNToGraph, PairsToDBN
    from .R3FUtils  import Base, BASE3, MODIFIED, AtomAtomDistance, BPS, BasePairReferenceDistance
    from .R3FUtils  import CoordsToMatrix, GAPS
    from .R3FReader import GuessFormat, ParsePDB, ParseCIF
    from .R3FWriter import WriteToPDB, WriteToCIF


def FileToEntry(filepath, fileformat = None):
    """Parse an input file into a PDB entry dictionary:
    Entry = {'IDLIST': [idstr1, idstr2, ...],
             'IDDICT': {idstr1: {'ATOMS' : ['N1','C2', ...],
                                 'N1': [X, Y, Z],
                                 ...}
                       ...}
            }
    """
    # Decide on the file format
    if not fileformat or fileformat.lower() not in {'pdb', 'cif', 'mmcif'}:
        fileformat = GuessFormat(filepath)
    else:
        fileformat = fileformat.lower().replace('mmcif','cif')

    # Parse the file
    try:
        if fileformat == 'pdb':
            Entry = ParsePDB(filepath)
        else:
            Entry = ParseCIF(filepath)
    except Exception as err:
        print("ERROR: failed to parse {}-format data from {}"\
              .format(fileformat.upper(),filepath))
        raise err

    return Entry


def EntryToFile(entry, filepath, fileformat = None):
    """ Print an input entry into a file (filepath)
        in PDB/mmCIF format (fileformat)
        By default, if the format is not specified by the user,
        it will be determined based on the file extension or,
        if it cannot be done, the format will be set to PDB."""

    #If None, get from the file extension
    if not fileformat:
        fileformat = filepath.split('.')[-1]

    # to lowercase
    fileformat = fileformat.lower()

    # if uncertain, set to pdb
    if fileformat not in ('pdb','cif','mmcif'):
        fileformat = 'pdb'

    try:
        if fileformat == 'pdb':
            WriteToPDB(entry, filepath)
        else:
            WriteToCIF(entry, filepath)
        return 0
    except Exception as err:
        print("ERROR: failed to print an entry in {} format to {}"\
              .format(fileformat.upper(), filepath))
        raise err


def CleanEntry(entry):
    """removes non-residues from the entry
    and returns the "clean" entry"""
    clean_entry = {}
    clean_entry['IDLIST'] = CleanIdList(entry)
    clean_entry['IDDICT'] = {idstr:entry['IDDICT'][idstr]
                             for idstr in clean_entry['IDLIST']}
    return clean_entry  


def CleanEntryToSequence(clean_entry, O3pP_dist = 2.0, O3p = "O3'", P = "P"):
    """Entry -> Sequence with 3'-ending residues in lowercase"""
    seq = ''

    L = len(clean_entry['IDLIST'])
    
    for i, idstr in enumerate(clean_entry['IDLIST']):

        base = Base(idstr)

        if base in MODIFIED:
            char = MODIFIED[base]
        else:
            char = 'N'

        # Next residue is i+1 unless we have the
        # last symbol in the sequence, than 0th
        # residue is considered the next to handle
        # circular molecules
        if i == L - 1:
            j = 0
        else:
            j = i + 1
        nxt = clean_entry['IDLIST'][j]
        
        if O3p in clean_entry['IDDICT'][idstr]\
           and P in clean_entry['IDDICT'][nxt]\
           and AtomAtomDistance(clean_entry['IDDICT'][idstr][O3p],
                                clean_entry['IDDICT'][nxt][P]) < O3pP_dist:
            seq += char
        else:
            seq += char.lower()

    return seq
    

def CleanEntryToBasePairs(clean_entry, thresh1 = 0.5, thresh2 = 2.0,
                          raw_filter = 15.0):
    """Annotate canonical GC, AU, and GU base pairs"""
    bps = []

    for i, idstr1 in enumerate(clean_entry['IDLIST']):
        for j, idstr2 in enumerate(clean_entry['IDLIST']):

            base1, base2 = Base(idstr1), Base(idstr2)

            if base1 in MODIFIED and base2 in MODIFIED:
                if MODIFIED[base1] + MODIFIED[base2] in BPS and\
                   all(BASE3[base1][k] in clean_entry['IDDICT'][idstr1] for k in range(3)) and\
                   all(BASE3[base2][k] in clean_entry['IDDICT'][idstr2] for k in range(3)) and\
                   AtomAtomDistance(clean_entry['IDDICT'][idstr1][BASE3[base1][0]],
                                    clean_entry['IDDICT'][idstr2][BASE3[base2][0]]) < raw_filter:
                    
                    bp = np.array([clean_entry['IDDICT'][idstr][k]
                                   for base,idstr in ((base1, idstr1),
                                                      (base2, idstr2))
                                   for k in BASE3[base]])

                    dist1, dist2 = BasePairReferenceDistance(bp,
                                                             MODIFIED[base1] + MODIFIED[base2])
                    if dist1 <= thresh1 and dist2 <= thresh2:
                        bps.append([dist1, min(i, j), max(i, j)])

    return sorted(bps)


def CleanEntryToSecStruct(clean_entry):

    seq     = CleanEntryToSequence(clean_entry)
    raw_bps = CleanEntryToBasePairs(clean_entry)
    bps     = []
    dbn     = ''

    paired  = set()

    for score,v,w in sorted(raw_bps):

        if v not in paired and w not in paired:
            bps.append((v, w))
            paired.add(v)
            paired.add(w)

    dbn = PairsToDBN(bps, len(seq))
    
    return seq, dbn, sorted(bps)


def EntryToStructure(entry):
    """ 
    Structure = {'IDLIST': [idstr1, idstr2, ...], <nucleic acid residues only>
                 'IDDICT': {idstr1: {'ATOMS' : ['N1','C2', ...],
                                                'N1': [X, Y, Z],
                                     ...}
                           ...},
                 'SEQ': string == unmodified sequence, 3'-ends in lowercase,
                 'DBN': string == dot-bracket notation (base pairs),
                 'BREAKS': list of binary 3'-ending/non-3'-ending residue values,
                 'MODIFIED': list of binary modified/standard residue values,
                 'BPS': list of base pairs, bp = [i, j, rmsd],
                 'GRAPH': {'ADJ': adjacency list (dictionary) {i: {j,...},},
                           'CC': connected components [{i,j,...},{...}]}
                 'CONNECTED': bool,
                 }
    """

    structure = CleanEntry(entry) # remove all "non-residues"

    seq, dbn, bps  = CleanEntryToSecStruct(structure)

    structure['SEQ'] = seq
    structure['DBN'] = dbn
    structure['BPS'] = bps
    
    structure['BREAKS']   = BreaksLine(structure['SEQ'])
    structure['MODIFIED'] = ModifiedLine(structure['IDLIST'])
    structure['GRAPH']    = SeqAndDBNToGraph(seq, dbn)

    if len(structure['GRAPH']['CC']) > 1:
        structure['CONNECTED'] = False
    else:
        structure['CONNECTED'] = True

    return structure

def FileToStructure(file):
    return EntryToStructure(FileToEntry(file))


def SplitStructures(struct):
    """structure -> list of connected structures"""
    if struct['CONNECTED']:
        return [struct,]

    structs = []

    for cc in struct['GRAPH']['CC']:

        newstruct = {}
        mapp = {} # old index : new index mapping

        newstruct['SEQ']      = ''
        newstruct['DBN']      = ''
        newstruct['IDLIST']   = []
        newstruct['BREAKS']   = []
        newstruct['MODIFIED'] = []

        cnt = 0
        for i in range(len(struct['SEQ'])):
            if i in cc:
                mapp[i] = cnt
                cnt += 1
                newstruct['SEQ'] += struct['SEQ'][i]
                newstruct['DBN'] += struct['DBN'][i]
                newstruct['IDLIST'].append(struct['IDLIST'][i])
                newstruct['BREAKS'].append(struct['BREAKS'][i])
                newstruct['MODIFIED'].append(struct['MODIFIED'][i])
                
        newstruct['IDDICT'] = {k:struct['IDDICT'][k] for k in newstruct['IDLIST']}
        newstruct['BPS'] = [(mapp[v],mapp[w]) for v,w in struct['BPS'] if v in cc]
        newstruct['CONNECTED'] = True
        newstruct['GRAPH'] = {'ADJ':{}, 'CC': []}
        newstruct['GRAPH']['CC'].append(set(range(len(cc))))

        for oldind in struct['GRAPH']['ADJ']:
            if oldind in cc:
                newstruct['GRAPH']['ADJ'][mapp[oldind]] = {mapp[w]:struct['GRAPH']['ADJ'][oldind][w]
                                                           for w in struct['GRAPH']['ADJ'][oldind]}
        structs.append(newstruct)
            
    return structs


def StructureToModel(struct, basesize = 3):
    """structure -> {'SEQ':seq, 'DBN':dbn, 'COORDS': nparray((3,3*len(seq))}"""

    seq    = struct['SEQ']
    dbn    = struct['DBN']
    coords = np.zeros((basesize*len(struct['SEQ']), 3))

    for i, idstr in enumerate(struct['IDLIST']):

        atomnames = BASE3[Base(idstr)]

        for mod in range(basesize):

            coords[basesize*i + mod] = struct['IDDICT'][idstr][atomnames[mod]]

    return {'SEQ': seq,
            'DBN': dbn,
            'COORDS': coords}


def ModelToMatrix(model):
    """model -> distance matrix"""
    return CoordsToMatrix(model['COORDS'])
    

def EntryToCoords(entry, basesize = 3):
    """entry -> nparray((3,3*len(entry['IDLIST']))"""

    coords = np.zeros((basesize*len(entry['IDLIST']), 3))

    for i, idstr in enumerate(entry['IDLIST']):

        atomnames = BASE3[Base(idstr)]

        for mod in range(basesize):
            coords[basesize*i + mod] = entry['IDDICT'][idstr][atomnames[mod]]

    return coords


def ParseRestraints(restraints, chunksymbol = "#"):
    """list of these:
       [filename, (target, template, dbn, ###line1, ###line2,..)]
       -->
       list of these: 
       {"IDLIST": [integers],
        "STRUCT": entry,
        "COORDS": N*3 np.array}
    """

    if not restraints:
        return None

    ###TODO - derive the missing lines
    ### for now - assume all are present
    ### if anything is missing - chage ToEntry to ToStructure

    ###TODO consistency checks is the lines are present (entry vs. lines)

    parsed_restraints = []

    for file in restraints:

        filename = file[0]
        target   = file[1]
        template = file[2]
        dbn      = file[3]
        chunks   = file[4:]
        
        entry = FileToEntry(file[0])

        for chunk in chunks:

            X, Y = [], []
            x, y = -1, -1

            for a,b,c in zip(template, target, chunk):
                if not a in GAPS:
                    x += 1
                if not b in GAPS:
                    y += 1
                if c == chunksymbol:
                    X.append(x)
                    Y.append(y)

            chunkentry = {'IDLIST': [],
                          'IDDICT': {}}
            for i in X:
                idstr = entry['IDLIST'][i]
                chunkentry['IDLIST'].append(idstr)
                chunkentry['IDDICT'][idstr] = entry['IDDICT'][idstr]

            parsed_restraints.append({'IDLIST': Y,
                                      'STRUCT': chunkentry,
                                      'COORDS': EntryToCoords(chunkentry)})
    
    return parsed_restraints
            

        

#####################################
if __name__ == "__main__":
    entry = FileToEntry("data/1ffk_0_kt7.cif")
    #entry = FileToEntry("data/3d2g.pdb")
    #entry = FileToEntry("data/2fmt.pdb")
    #entry = FileToEntry("data/4v8o.cif")
    struct = EntryToStructure(entry)
    structs = SplitStructures(struct)
    for struct in structs:
        model = StructureToModel(struct)
        print(model['SEQ'])
        print(model['DBN'])
        print(model['COORDS'])
        print(ModelToMatrix(model))
#####################################
