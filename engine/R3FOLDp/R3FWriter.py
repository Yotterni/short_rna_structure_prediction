try:
    from R3FUtils import IntToChain
except:
    from .R3FUtils import IntToChain

TYPE_SYMBOLS_SRC = '''

AG	AG
AS	AS
AU	AU
BA	BA
BR	BR
BR3	BR
C	C
C01	C
C02	C
C03	C
C04	C
C05	C
C06	C
C07	C
C08	C
C09	C
C0B	C
C1	C
C1'	C
C1'1	C
C10	C
C11	C
C118	C
C12	C
C120	C
C121	C
C13	C
C14	C
C140	C
C141	C
C143	C
C144	C
C15	C
C16	C
C17	C
C18	C
C19	C
C1A	C
C1B	C
C1C	C
C1D	C
C1E	C
C1F	C
C1G	C
C1H	C
C1I	C
C1J	C
C1K	C
C1L	C
C1M	C
C1N	C
C1O	C
C1P	C
C1Q	C
C1R	C
C1S	C
C1T	C
C1U	C
C1V	C
C1W	C
C1X	C
C1Y	C
C2	C
C2'	C
C2'1	C
C20	C
C21	C
C22	C
C23	C
C24	C
C25	C
C26	C
C27	C
C28	C
C29	C
C2A	C
C2B	C
C2C	C
C2D	C
C2E	C
C2G	C
C2M	C
C2N	C
C2P	C
C2R	C
C3	C
C3'	C
C3'1	C
C30	C
C31	C
C32	C
C33	C
C34	C
C35	C
C36	C
C37	C
C38	C
C39	C
C3A	C
C3B	C
C3C	C
C3D	C
C3E	C
C3G	C
C3N	C
C3P	C
C3R	C
C3U	C
C4	C
C4'	C
C4'1	C
C40	C
C41	C
C42	C
C43	C
C44	C
C45	C
C46	C
C47	C
C48	C
C49	C
C4A	C
C4B	C
C4C	C
C4D	C
C4E	C
C4N	C
C4R	C
C5	C
C5'	C
C5'1	C
C50	C
C51	C
C52	C
C53	C
C54	C
C55	C
C56	C
C57	C
C58	C
C59	C
C5A	C
C5B	C
C5C	C
C5D	C
C5E	C
C5M	C
C5N	C
C5R	C
C6	C
C6'	C
C60	C
C61	C
C62	C
C63	C
C64	C
C65	C
C67	C
C68	C
C69	C
C6A	C
C6B	C
C6C	C
C6M	C
C6N	C
C7	C
C7'	C
C70	C
C71	C
C72	C
C73	C
C74	C
C75	C
C76	C
C77	C
C79	C
C7A	C
C7B	C
C7M	C
C7N	C
C7X	C
C8	C
C80	C
C81	C
C82	C
C83	C
C84	C
C85	C
C86	C
C87	C
C8A	C
C8B	C
C8C	C
C8M	C
C9	C
C91	C
C92	C
C93	C
C94	C
C9A	C
C9B	C
CA	CA
CA'	C
CA1	C
CA2	C
CA3	C
CA4	C
CA5	C
CA6	C
CA7	C
CA8	C
CA9	C
CAA	C
CAB	C
CAC	C
CAD	C
CAE	C
CAF	C
CAG	C
CAH	C
CAI	C
CAJ	C
CAK	C
CAL	C
CAM	C
CAN	C
CAO	C
CAP	C
CAQ	C
CAR	C
CAS	C
CAT	C
CAU	C
CAV	C
CAW	C
CAX	C
CAY	C
CAZ	C
CB	C
CB'	C
CB1	C
CB2	C
CB3	C
CB4	C
CB5	C
CB6	C
CBA	C
CBB	C
CBC	C
CBD	C
CBE	C
CBF	C
CBG	C
CBH	C
CBI	C
CBK	C
CBM	C
CBN	C
CBO	C
CBP	C
CBQ	C
CBR	C
CBS	C
CBT	C
CBU	C
CBV	C
CBW	C
CBX	C
CBY	C
CBZ	C
CC	C
CC1	C
CC2	C
CC3	C
CC4	C
CC5	C
CC6	C
CD	C
CD1	C
CD2	C
CE	C
CE1	C
CE2	C
CE3	C
CG	C
CG1	C
CG2	C
CH2	C
CH3	C
CHA	C
CHB	C
CHC	C
CHD	C
CI3	C
CL	CL
CL1	CL
CL7	CL
CM	C
CM'	C
CM1	C
CM2	C
CM4	C
CM5	C
CM6	C
CM7	C
CN	C
CN1	C
CN7	C
CO	CO
CS	CS
CV'	C
CW'	C
CX'	C
CX3	C
CX4	C
CX5	C
CX6	C
CXD	C
CXE	C
CXF	C
CXG	C
CXN	C
CXO	C
CXP	C
CXQ	C
CY'	C
CZ	C
CZ'	C
CZ2	C
CZ3	C
F	F
F01	F
F1	F
F1'	F
F2	F
F2'	F
F20	F
F21	F
F22	F
F3	F
F3'	F
F4	F
F5	F
F61	F
F7	F
F8	F
FAI	F
FE	FE
H	H
H1	H
H1'	H
H1'1	H
H10	H
H101	H
H102	H
H103	H
H11	H
H111	H
H112	H
H12	H
H121	H
H122	H
H13	H
H131	H
H132	H
H14	H
H15	H
H16	H
H161	H
H162	H
H163	H
H17	H
H171	H
H172	H
H18	H
H181	H
H182	H
H183	H
H19	H
H1A	H
H1B	H
H1D	H
H2	H
H2'	H
H2''	H
H2'1	H
H20	H
H21	H
H22	H
H23	H
H24	H
H25	H
H26	H
H261	H
H262	H
H263	H
H27	H
H271	H
H272	H
H28	H
H281	H
H282	H
H283	H
H29	H
H2A	H
H2B	H
H2C	H
H2D	H
H2N	H
H2O1	H
H3	H
H3'	H
H3'1	H
H3'2	H
H31	H
H32	H
H33	H
H34	H
H35	H
H36	H
H361	H
H362	H
H363	H
H37	H
H371	H
H372	H
H38	H
H381	H
H382	H
H39	H
H3A	H
H3B	H
H3D	H
H3U1	H
H3U2	H
H3U3	H
H4	H
H4'	H
H4'1	H
H40	H
H41	H
H42	H
H45	H
H461	H
H462	H
H463	H
H471	H
H472	H
H481	H
H482	H
H4A	H
H4B	H
H4D	H
H4N	H
H5	H
H5'	H
H5''	H
H5'1	H
H5'2	H
H5'A	H
H51	H
H51A	H
H51N	H
H52	H
H52A	H
H52N	H
H5A1	H
H5A2	H
H5B1	H
H5B2	H
H5M	H
H5N	H
H6	H
H61	H
H61A	H
H62	H
H62A	H
H6L	H
H6M	H
H6N	H
H7	H
H71	H
H71N	H
H72	H
H72N	H
H73	H
H8	H
H81	H
H82	H
H83	H
H8A	H
H8C	H
H9	H
H91	H
H92	H
H93	H
H9C1	H
H9C2	H
HA	H
HA2	H
HA3	H
HAA	H
HAC	H
HAE	H
HAF	H
HAG	H
HAH	H
HAI	H
HAJ	H
HAK	H
HAL	H
HAM	H
HAN	H
HAO	H
HAP	H
HAR	H
HAS	H
HAT	H
HAU	H
HB	H
HB1	H
HB2	H
HB3	H
HD1	H
HD11	H
HD12	H
HD13	H
HD2	H
HD21	H
HD22	H
HD23	H
HD3	H
HE	H
HE1	H
HE2	H
HE21	H
HE22	H
HE3	H
HG	H
HG1	H
HG11	H
HG12	H
HG13	H
HG2	H
HG21	H
HG22	H
HG23	H
HG3	H
HH	H
HH11	H
HH12	H
HH2	H
HH21	H
HH22	H
HM'1	H
HM'2	H
HM'3	H
HM11	H
HM12	H
HM13	H
HM21	H
HM22	H
HM23	H
HM51	H
HM52	H
HM53	H
HM71	H
HM72	H
HM73	H
HN0	H
HN1	H
HN10	H
HN11	H
HN12	H
HN13	H
HN2	H
HN21	H
HN22	H
HN23	H
HN3	H
HN31	H
HN32	H
HN33	H
HN4	H
HN41	H
HN42	H
HN43	H
HN5	H
HN51	H
HN52	H
HN53	H
HN6	H
HN61	H
HN62	H
HN63	H
HN7	H
HN91	H
HO1	H
HO2	H
HO2'	H
HO2A	H
HO2B	H
HO2N	H
HO3	H
HO3'	H
HO3A	H
HO3B	H
HO3N	H
HO4	H
HO5	H
HO5'	H
HO6	H
HZ	H
HZ1	H
HZ2	H
HZ3	H
I	I
I5	I
IR	IR
K	K
LU	LU
MG	MG
MN	MN
N	N
N01	N
N02	N
N03	N
N04	N
N05	N
N06	N
N07	N
N08	N
N09	N
N0A	N
N1	N
N1'	N
N10	N
N11	N
N12	N
N13	N
N139	N
N14	N
N15	N
N16	N
N17	N
N19	N
N1A	N
N1B	N
N1C	N
N1N	N
N2	N
N2'	N
N20	N
N21	N
N22	N
N23	N
N24	N
N25	N
N26	N
N27	N
N28	N
N29	N
N2A	N
N2B	N
N3	N
N3'	N
N31	N
N32	N
N33	N
N34	N
N35	N
N37	N
N39	N
N3A	N
N3B	N
N3C	N
N4	N
N4'	N
N40	N
N41	N
N45	N
N49	N
N4A	N
N5	N
N50	N
N51	N
N52	N
N53	N
N58	N
N59	N
N6	N
N61	N
N62	N
N63	N
N64	N
N6A	N
N6C	N
N7	N
N71	N
N72	N
N7A	N
N7B	N
N7C	N
N7N	N
N8	N
N82	N
N9	N
N91	N
N9A	N
N9B	N
N9C	N
NA	NA
NA2	N
NA7	N
NAA	N
NAB	N
NAC	N
NAE	N
NAG	N
NAH	N
NAK	N
NAL	N
NAM	N
NAR	N
NAS	N
NAT	N
NAW	N
NB	N
NB1	N
NB4	N
NBD	N
NBE	N
NBF	N
NBG	N
NBH	N
NBI	N
NBJ	N
NBK	N
NBL	N
NBN	N
NBP	N
NBQ	N
NBV	N
NC	N
NC1	N
NC4	N
NC6	N
NCA	N
NCB	N
NCC	N
NCD	N
ND	N
ND1	N
ND2	N
NE	N
NE1	N
NE2	N
NF1	N
NH1	N
NH2	N
NI	NI
NXT	N
NXU	N
NXV	N
NXW	N
NZ	N
NZ'	N
O	O
O01	O
O02	O
O03	O
O04	O
O05	O
O06	O
O07	O
O08	O
O09	O
O1	O
O1'	O
O10	O
O11	O
O12	O
O13	O
O14	O
O142	O
O15	O
O16	O
O17	O
O18	O
O19	O
O1A	O
O1B	O
O1C	O
O1D	O
O1G	O
O1N	O
O1S	O
O1V	O
O1X	O
O2	O
O2'	O
O2'1	O
O20	O
O21	O
O22	O
O23	O
O24	O
O25	O
O27	O
O28	O
O29	O
O2A	O
O2B	O
O2C	O
O2D	O
O2E	O
O2G	O
O2N	O
O2R	O
O2S	O
O2V	O
O2X	O
O3	O
O3'	O
O3'1	O
O30	O
O31	O
O32	O
O33	O
O34	O
O35	O
O36	O
O39	O
O3A	O
O3B	O
O3C	O
O3D	O
O3E	O
O3G	O
O3P	O
O3R	O
O3S	O
O3X	O
O4	O
O4'	O
O4'1	O
O40	O
O41	O
O42	O
O43	O
O44	O
O45	O
O48	O
O4A	O
O4B	O
O4D	O
O4E	O
O4P	O
O4R	O
O5	O
O5'	O
O5'1	O
O50	O
O51	O
O52	O
O53	O
O54	O
O56	O
O58	O
O5A	O
O5B	O
O5D	O
O5E	O
O5P	O
O5R	O
O6	O
O61	O
O62	O
O63	O
O66	O
O6A	O
O6B	O
O6P	O
O6R	O
O7	O
O72	O
O74	O
O78	O
O7A	O
O7N	O
O7R	O
O8	O
O80	O
O83	O
O85	O
O86	O
O8A	O
O8R	O
O9	O
OA1	O
OA4	O
OA5	O
OA6	O
OA8	O
OAB	O
OAC	O
OAD	O
OAE	O
OAF	O
OAG	O
OAH	O
OAI	O
OAS	O
OAT	O
OAU	O
OAV	O
OAX	O
OAY	O
OAZ	O
OB1	O
OB2	O
OB3	O
OB6	O
OBA	O
OBB	O
OBD	O
OBE	O
OBI	O
OBJ	O
OBL	O
OBR	O
OBU	O
OC	O
OC1	O
OC2	O
OC3	O
OCA	O
OCB	O
OD1	O
OD2	O
ODA	O
ODB	O
OE1	O
OE2	O
OG	O
OG1	O
OG2	O
OH	O
OH2	O
OH3	O
OH4	O
OH5	O
OH6	O
OH7	O
OM5	O
OM7	O
ON1	O
ON2	O
OO	O
OP	O
OP1	O
OP11	O
OP2	O
OP21	O
OP3	O
OP4	O
OP5	O
OP6	O
OS	OS
OW'	O
OX'	O
OXT	O
OY'	O
P	P
P'	P
P01	P
P02	P
P03	P
P1	P
P11	P
P2	P
P22	P
P26	P
P3	P
P30	P
PA	P
PAX	P
PAZ	P
PB	P
PBK	P
PBR	P
PBS	P
PC	P
PD	P
PG	P
PN	P
RH	RH
S	S
S1	S
S10	S
S15	S
S2	S
S2'	S
S23	S
S26	S
S2P	S
S3P	S
S4	S
S4'	S
S6	S
SD	S
SE	SE
SE1	SE
SE2	SE
SE2'	SE
SE4	SE
SG	S
SM	SM
SP1	S
SP2	S
SR	SR
TB	TB
TL	TL
U	U
UNK	X
V	V
ZN	ZN
'''

TYPE_SYMBOLS = {x.split()[0]:x.split()[1] for x in TYPE_SYMBOLS_SRC.strip().split('\n')}


def Atom(idstr):
    '''"1.A.G.22.B" -> Atom dict:
       atom = {"pdbx_PDB_model_num" : '1',
               "auth_asym_id"       : 'A',
               "auth_comp_id"       : 'G',
               "auth_seq_id"        : '22',
               "pdbx_PDB_ins_code"  : 'B',
              }
    '''

    splt = idstr.split('.')
    assert len(splt) == 5, "Incorrect DSSR-like id string: {}".format(idstr)

    atom = {"pdbx_PDB_model_num": splt[0],
            "auth_asym_id"      : splt[1],
            "auth_comp_id"      : splt[2],
            "auth_seq_id"       : splt[3],
            "pdbx_PDB_ins_code" : splt[4],
            }

    return atom


def WriteAtomToPDB(atom):
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    -------------------------------------------------------------------------------------
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
    """

    res = "ATOM  "
    atomid = str(atom['id'] % 10**5)
    res += ' '*(5 - len(atomid)) + atomid
    res += ' '

    atom['type_symbol'] = TYPE_SYMBOLS[atom['auth_atom_id']]

    if len(atom["auth_atom_id"]) == 4:
        res += atom["auth_atom_id"]
    else:
        half1 = atom["type_symbol"]
        half2 = atom["auth_atom_id"][len(atom["type_symbol"]):]
        res += ' '*(2 - len(half1)) + half1
        res += half2 + ' '*(2-len(half2))

    res += ' ' #label_alt_id

    res += ' '*(3 - len(atom["auth_comp_id"])) + atom["auth_comp_id"]

    res += ' '*(2 - len(atom["auth_asym_id"])) + atom["auth_asym_id"]

    residue_id = str(int(atom["auth_seq_id"]) % 10**4)
    res += ' '*(4 - len(residue_id)) + residue_id

    res += ' '*(1 - len(atom["pdbx_PDB_ins_code"])) + atom["pdbx_PDB_ins_code"]

    res += ' ' * 3
    
    for Cartn_coord in ("Cartn_x", "Cartn_y", "Cartn_z"):
        coord = "{:.3f}".format(round(atom[Cartn_coord], 3))
        res += ' '*(8 - len(coord)) + coord

    occup = "{:.2f}".format(round(atom["occupancy"], 2))
    res += ' '*(6 - len(occup)) + occup

    bfactor = "{:.2f}".format(round(atom["B_iso_or_equiv"], 2))
    res += ' '*(6 - len(bfactor)) + bfactor

    res += ' ' * 10

    res += ' '*(2 - len(atom["type_symbol"])) + atom["type_symbol"]

    res += ' ' * 2 #pdbx_formal_charge

    return res


def WriteTER(atom):

    res = ""
    res += "TER" + " "*3

    atomid = str(atom['id'] % 10**5)
    res += ' '*(5 - len(atomid)) + atomid

    res += ' ' * 6

    res += ' '*(3 - len(atom["auth_comp_id"])) + atom["auth_comp_id"]

    res += ' '*(2 - len(atom["auth_asym_id"])) + atom["auth_asym_id"]

    residue_id = str(int(atom["auth_seq_id"]) % 10**4)
    res += ' '*(4 - len(residue_id)) + residue_id

    res += ' '*(1 - len(atom["pdbx_PDB_ins_code"])) + atom["pdbx_PDB_ins_code"]

    return res
    

def WriteToPDB(entry, file):

    if not entry['IDLIST']:
        raise Exception("ERROR: entry IDLIST field is empty (WriteToPDB)")

    # get model identifiers and sort them
    models = set()
    for idstr in entry['IDLIST']:
        try:
            models.add(int(Atom(idstr)["pdbx_PDB_model_num"]))
        except:
            pass # ignore non-int models, set them to 1 later
    if not models: # if no integer models, consider everything model 1
        models.add(1) 
    models = [str(m) for m in sorted(models)]

    # Handle chain identifiers longer than 2-letter
    # By renaming all chains from A to zz
    newchains = {}
    if any(len(Atom(idstr)["auth_asym_id"]) > 2
           for idstr in entry['IDLIST']):
        chaincount = 0
        for idstr in entry['IDLIST']:
            chain = Atom(idstr)["auth_asym_id"]
            if chain not in newchains:
                newchains[chain] = IntToChain(chaincount)
                chaincount += 1        

    with open(file, 'w') as outp:

        # if renaming chains
        if newchains:
            outp.write("REMARK 250 RENAMED CHAINS:" + '\n')
            for chain in sorted(newchains.keys()):
                outp.write("REMARK 250 {} -> {}".format(chain, newchains[chain]) + '\n')
        
        for model in models:

            atomcount = 1
            prevatom  = None
            outp.write("MODEL        {} \n".format(model))

            for idstr in entry['IDLIST']:
                atom = Atom(idstr)

                # default occupancy and B-factor
                atom["occupancy"]      = 1.00 
                atom["B_iso_or_equiv"] = 0.00

                # if renaming chains
                if newchains:
                    atom["auth_asym_id"] = newchains[atom["auth_asym_id"]]
                
                # if uncertain model -> set to 1
                if not atom["pdbx_PDB_model_num"].isdigit():
                    atom["pdbx_PDB_model_num"] = '1'
                # skip residues from other models
                if atom["pdbx_PDB_model_num"] != model:
                    continue

                # Print TER line if new chain
                if atomcount > 1:
                    if atom['auth_asym_id'] != prevatom['auth_asym_id']:
                        outp.write(WriteTER(prevatom)+'\n')
                
                for atom_id in entry['IDDICT'][idstr]['ATOMS']:
                            
                    atom['id'] = atomcount
                    atomcount += 1

                    atom['auth_atom_id'] = atom_id
                    atom['Cartn_x'] = entry['IDDICT'][idstr][atom_id][0]
                    atom['Cartn_y'] = entry['IDDICT'][idstr][atom_id][1]
                    atom['Cartn_z'] = entry['IDDICT'][idstr][atom_id][2]

                    outp.write(WriteAtomToPDB(atom)+'\n')

                prevatom = {k:v for k,v in atom.items()}
      
            # Print TER line after the last atom in the model
            outp.write(WriteTER(atom)+'\n')

            outp.write("ENDMDL\n")
        outp.write("END\n")          



def WriteToCIF(entry, file):

    if not entry['IDLIST']:
        raise Exception("ERROR: entry IDLIST field is empty (WriteToCIF)")

    title = ["group_PDB",
             "id",
             "type_symbol",
             "label_atom_id",
             "label_alt_id",
             "label_comp_id",
             "label_asym_id",
             "label_entity_id",
             "label_seq_id",
             "pdbx_PDB_ins_code",
             "Cartn_x",
             "Cartn_y",
             "Cartn_z",
             "occupancy",
             "B_iso_or_equiv",
             "pdbx_formal_charge",
             "auth_seq_id",
             "auth_comp_id",
             "auth_asym_id",
             "auth_atom_id",
             "pdbx_PDB_model_num",
             ]

    atomcount    = 1
    residuecount = 1

    entities = {} # derive entity_id as index of the (model, chain) pair

    with open(file, 'w') as outp:

        outp.write("#\nloop_\n")

        for col in title:
            outp.write("_atom_site."+col+'\n')

        for idstr in entry['IDLIST']:
            atom = Atom(idstr)

            atom["group_PDB"] = "ATOM"
            atom["label_seq_id"] = residuecount
            residuecount += 1

            # default occupancy and B-factor
            atom["occupancy"]      = 1.00 
            atom["B_iso_or_equiv"] = 0.00

            # if uncertain model -> set to 1
            if not atom["pdbx_PDB_model_num"].isdigit():
                atom["pdbx_PDB_model_num"] = '1'

            if (atom["pdbx_PDB_model_num"],
                atom["auth_asym_id"]) not in entities:
                entity_id = len(entities) + 1
                entities[(atom["pdbx_PDB_model_num"],
                          atom["auth_asym_id"])] = entity_id

            atom["label_entity_id"] = entities[(atom["pdbx_PDB_model_num"],
                                                atom["auth_asym_id"])]
            atom["label_alt_id"]  = '.'
            atom["label_comp_id"] = atom["auth_comp_id"]
            atom["label_asym_id"] = atom["auth_asym_id"]

            # avoid empty strings
            for col in title:
                if col not in atom or not str(atom[col]).strip():
                    atom[col] = '?'
                
            for atom_id in entry['IDDICT'][idstr]['ATOMS']:
                            
                atom['id'] = atomcount
                atomcount += 1

                atom['auth_atom_id']  = atom_id
                atom["label_atom_id"] = atom['auth_atom_id']
                atom['type_symbol']   = TYPE_SYMBOLS[atom['auth_atom_id']]
                atom['Cartn_x'] = entry['IDDICT'][idstr][atom_id][0]
                atom['Cartn_y'] = entry['IDDICT'][idstr][atom_id][1]
                atom['Cartn_z'] = entry['IDDICT'][idstr][atom_id][2]

                outp.write(' '.join([str(atom[col]) for col in title]) + '\n')

        outp.write("#\n")
