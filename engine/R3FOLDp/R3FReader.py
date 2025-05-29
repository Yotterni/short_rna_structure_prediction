
def DSSRidstr(atom):
    '''1.A.G.22.B - DSSR-like id string
       [model.chain.res.resnum.inscode]'''
    return '.'.join([str(x) for x in [atom["pdbx_PDB_model_num"],
                                      atom["auth_asym_id"],
                                      atom["auth_comp_id"],
                                      atom["auth_seq_id"],
                                      atom["pdbx_PDB_ins_code"]]])

def GuessFormat(filepath):
    """Determine the format of a file -> 'cif' or 'pdb' """

    pdb = 0
    cif = 0

    with open(filepath) as file:
        for line in file:
            if line.startswith('#'):
                cif += 1
            elif line.startswith('_atom.site'):
                cif += 1
            elif line.startswith('loop_'):
                cif += 1
            elif line.startswith('END'):
                pdb += 1
            elif line.startswith('MODEL'):
                pdb += 1
            elif line.startswith('TER'):
                pdb += 1
            elif line.startswith('CONECT'):
                pdb += 1   
                
    return ['pdb','cif'][int(cif > pdb)]

def ParseAtomPDB(line):
    '''https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html'''
    
    atom = {}
    atom["group_PDB"]          = line[:6].strip()
    atom["id"]                 = int(line[6:11])
    atom["auth_atom_id"]       = line[12:16].strip()
    atom["label_alt_id"]       = line[16].strip()
    atom["auth_comp_id"]       = line[17:20].strip()
    atom["auth_asym_id"]       = line[20:22].strip()
    atom["auth_seq_id"]        = int(line[22:26])
    atom["pdbx_PDB_ins_code"]  = line[26].strip()
    atom["Cartn_x"]            = float(line[30:38])
    atom["Cartn_y"]            = float(line[38:46])
    atom["Cartn_z"]            = float(line[46:54])
    atom["occupancy"]          = float(line[54:60]) if line[54:60].strip() else float("nan")
    atom["B_iso_or_equiv"]     = float(line[60:66]) if line[60:66].strip() else float("nan")
    atom["type_symbol"]        = line[76:78].strip()
    atom["pdbx_formal_charge"] = line[78:80].strip()

    # Convert SPDBV asterisks into prime symbols
    atom["auth_atom_id"] = atom["auth_atom_id"].replace('*',"'")
    # Convert SPDBV phosphate oxygen format O1P -> OP1, O2P -> OP2
    atom["auth_atom_id"] = atom["auth_atom_id"].replace('O1P','OP1').replace('O2P','OP2')

    if atom["pdbx_PDB_ins_code"] == '?':
        atom["pdbx_PDB_ins_code"] = ''

    if atom["pdbx_formal_charge"] == '?':
        atom["pdbx_formal_charge"] = ''

    if atom["label_alt_id"] == '.':
        atom["label_alt_id"] = ''

    return atom

def ParsePDB(filepath):
    '''https://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html#ATOMP'''

    # By default, the model is set to 1
    model  = "1"

    title1 = ["group_PDB","id","auth_atom_id",
             "label_alt_id","auth_comp_id","auth_asym_id",
             "auth_seq_id","pdbx_PDB_ins_code","Cartn_x",
             "Cartn_y","Cartn_z","occupancy",
             "B_iso_or_equiv","type_symbol","pdbx_formal_charge"]
    title2 = ["label_asym_id","label_atom_id","label_comp_id",
              "label_entity_id","label_seq_id","pdbx_PDB_model_num"]

    titlei = {x:i for i,x in enumerate(title1+title2)}

    residuelst, residuedct = [], {}

    with open(filepath) as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):

                atom = ParseAtomPDB(line)

                atom["pdbx_PDB_model_num"] = model
                atom["label_asym_id"] = atom["auth_asym_id"]
                atom["label_atom_id"] = atom["auth_atom_id"]
                atom["label_comp_id"] = atom["auth_comp_id"]
                atom["label_seq_id"]  = atom["auth_seq_id"]
                atom['DSSR']          = DSSRidstr(atom)

                # new residue found
                if atom['DSSR'] not in residuedct:
                    residuelst.append(atom['DSSR'])
                    residuedct[atom['DSSR']] = {'ATOMS': []}
                    
                residue = residuedct[atom['DSSR']]
                # new atom found (otherwise we have a new altloc for an already found atom)
                if atom['auth_atom_id'] not in residue:
                    residue['ATOMS'].append(atom['auth_atom_id'])
                # consider only the last altloc version of the atom and ignore all the previous ones
                residue[atom['auth_atom_id']] = [atom[t] for t in ["Cartn_x", "Cartn_y", "Cartn_z"]]
                
            elif line.startswith('MODEL'):
                model = line.strip().split()[-1]

    Entry = {'IDLIST': residuelst,
             'IDDICT': residuedct}

    return Entry


def ParseAtomCIF(line, title):

    linesplit = line.strip().split()
    atom = {title[i]:linesplit[i] for i in range(len(title))}

    for frag in ("atom", "comp", "asym", "seq"):

        auth  =  "auth_{}_id".format(frag)
        label = "label_{}_id".format(frag)
        
        if auth not in atom and label in atom:
            atom[auth] = atom[label]
        elif label not in atom and auth in atom:
            atom[label] = atom[auth]

    # Convert all integers from strings
    for int_token in ("id", "auth_seq_id"):
        atom[int_token] = int(atom[int_token]) if int_token in atom else float("nan")

    # By default, the model is set to 1
    atom["pdbx_PDB_model_num"] = int(atom["pdbx_PDB_model_num"])\
                                 if "pdbx_PDB_model_num" in atom and\
                                 atom["pdbx_PDB_model_num"].isdigit() else 1

    # Convert all floats from strings
    for float_token in ("Cartn_x", "Cartn_y", "Cartn_z","occupancy","B_iso_or_equiv"):
        atom[float_token] = float(atom[float_token]) if float_token in atom else float("nan")
    
    if      "auth_atom_id"  not in atom: atom["auth_atom_id"]       = ''
    if     "label_atom_id"  not in atom: atom["label_atom_id"]      = ''
    if      "label_alt_id"  not in atom: atom["label_alt_id"]       = ''
    if      "auth_comp_id"  not in atom: atom["auth_comp_id"]       = ''
    if     "label_comp_id"  not in atom: atom["label_comp_id"]      = ''
    if      "auth_asym_id"  not in atom: atom["auth_asym_id"]       = ''
    if "pdbx_PDB_ins_code"  not in atom: atom["pdbx_PDB_ins_code"]  = ''
    if "pdbx_formal_charge" not in atom: atom["pdbx_formal_charge"] = ''

    if atom["pdbx_PDB_ins_code"] == '?':
        atom["pdbx_PDB_ins_code"] = ''

    if atom["pdbx_formal_charge"] == '?':
        atom["pdbx_formal_charge"] = ''

    if atom["label_alt_id"] == '.':
        atom["label_alt_id"] = ''

    # Strip double quotes for cases like this: "O3'"
    atom["auth_atom_id"] = atom["auth_atom_id"].strip('"') 
    atom["label_atom_id"] = atom["label_atom_id"].strip('"')

    # Convert SPDBV asterisks into prime symbols
    atom["auth_atom_id"]  = atom["auth_atom_id"].replace('*',"'")
    atom["label_atom_id"] = atom["label_atom_id"].replace('*',"'")
    # Convert SPDBV phosphate oxygen format O1P -> OP1, O2P -> OP2
    atom["auth_atom_id"]  = atom["auth_atom_id"].replace('O1P','OP1').replace('O2P','OP2')
    atom["label_atom_id"] = atom["label_atom_id"].replace('O1P','OP1').replace('O2P','OP2')

    return atom


def ParseCIF(filepath):
    
    title  = []
    titlei = None

    residuelst, residuedct = [], {}

    with open(filepath) as file:
        for line in file:
            if line.startswith("_atom_site."):
                title.append(line.strip().split('.')[-1])

            elif line.startswith('ATOM') or line.startswith('HETATM'):

                atom = ParseAtomCIF(line, title)

                if not titlei:

                    title2 = [k for k in atom if k not in title]
                    finaltitle  = title + title2
                    titlei = {x:i for i,x in enumerate(finaltitle)}

                atom['DSSR'] = DSSRidstr(atom)

                # new residue found
                if atom['DSSR'] not in residuedct:
                    residuelst.append(atom['DSSR'])
                    residuedct[atom['DSSR']] = {'ATOMS': []}

                residue = residuedct[atom['DSSR']]
                # new atom found (otherwise we have a new altloc for an already found atom)
                if atom['auth_atom_id'] not in residue:
                    residue['ATOMS'].append(atom['auth_atom_id'])
                # consider only the last altloc version of the atom and ignore all the previous ones
                residue[atom['auth_atom_id']] = [atom[t] for t in ["Cartn_x", "Cartn_y", "Cartn_z"]]
                
            elif line.startswith('MODEL'):
                model = line.strip().split()[-1]

    Entry = {'IDLIST': residuelst,
             'IDDICT': residuedct}

    return Entry


