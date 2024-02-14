## BASIC LIBS
import os
from os import path as p
import pandas as pd


def right_aligned(pdbList, textList):
    lastSpaceIndex = len(pdbList) - 1 - pdbList[::-1].index(" ")
    # Replace spaces with letters starting from the end
    for i, letter in enumerate(reversed(textList)):
        pdbList[lastSpaceIndex - i] = letter
    return pdbList

def left_aligned(pdbList, textList):
    for i, letter in enumerate(textList):
        pdbList[i] = letter
    return pdbList


########################## write pdb dataframe to pdb file
def df2Pdb(df, outFile,
           chain=True):
    with open(outFile,"w") as f:
        for _, row in df.iterrows():

            pdbList = [" " for _ in range(80)]
            pdbList[0:6]    = left_aligned(pdbList[0:6],list(row['ATOM']))
            pdbList[6:11]   = right_aligned(pdbList[6:11],list(str(row['ATOM_ID'])))
            pdbList[12:16]  = right_aligned(pdbList[12:16], list(row['ATOM_NAME']))
            pdbList[17:20]  = right_aligned(pdbList[17:20], list(row['RES_NAME']))
            try: 
                pdbList[21]  = row['CHAIN_ID']
            except:
                pdbList[21] = " "
            pdbList[22:26]  = right_aligned(pdbList[22:26], list(str(row['RES_ID'])))
            pdbList[30:38]  = right_aligned(pdbList[30:38], list(f"{row['X']:>8.3f}"))
            pdbList[38:46]  = right_aligned(pdbList[38:46], list(f"{row['Y']:>8.3f}"))
            pdbList[46:54]  = right_aligned(pdbList[46:54], list(f"{row['Z']:>8.3f}"))
            pdbList[54:60]  = right_aligned(pdbList[54:60], list(f"{row['OCCUPANCY']:>6.2f}"))
            try:
                pdbList[60:66]  = right_aligned(pdbList[60:66], list(f"{row['BETAFACTOR']:>6.2f}"))
            except:
                pdbList[60:66] = right_aligned(pdbList[60:66],list("1.00"))
            try:
                pdbList[76:78] = right_aligned(pdbList[76:78], list(row['ELEMENT']))
            except:
                element = row["ATOM_NAME"][0]
                pdbList[76:78] = right_aligned(pdbList[76:78], list(element))
            pdbLine = "".join(pdbList)
            f.write(f"{pdbLine}\n")
        f.write("TER")
############################### read pdb file as pdb dataframe
def pdb2df(protPdb):
    columns = ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']
    data = []
    with open(protPdb, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_type = line[0:6].strip()
                try:
                    atom_id = int(line[6:11].strip())
                except:
                    atom_id = str(line[6:11].strip())

                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                if chain_id == '':
                    chain_id = None
                try:
                    res_id = int(line[22:26].strip())
                except:
                    res_id = str(line[22:26].strip())

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                try:
                    temp_factor = float(line[60:66].strip())
                except:
                    temp_factor = " "*6
                element = line[76:78].strip()

                data.append([atom_type, atom_id, atom_name, res_name, chain_id, res_id, x, y, z, occupancy, temp_factor, element])

    return pd.DataFrame(data, columns=columns)
############################### read list of pdb files as dataframes, concatonate, then write back to pdb file
def mergePdbs(pdbList,outFile):
    dfList=[]
    for pdbFile in pdbList:
        df = pdb2df(pdbFile)
        dfList.append(df)
    mergedDf = pd.concat(dfList,axis=0)
    df2Pdb(df=mergedDf, outFile=outFile, chain=True)

############################### apply a a bunch of fixes to a pdb dataframe
def fix_atom_names(df): 
    # deal with unwanted apostrophies (prime)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("'", "")
    # deal with numbers at the beginnig of atom names
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].replace(r'^(\d+)(.+)$', r'\2\1', regex=True)
    # deal with "A" at the start of atom name
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].apply(lambda x: x.lstrip('A') if x.startswith('A') else x)

    ## ensure unique names
    count_series = df.groupby('ATOM_NAME').cumcount()
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'] + "_" +count_series.astype(str)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_0", "")
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_", "")

    return df 

############################# reads a pdbqt file to pandas dataframe
def pdbqt2df(pdbqtFile):
    # remove ROOT/BRANCH    
    pdbqtColumns    =   ["ATOM","ATOM_ID", "ATOM_NAME", "RES_NAME",
                    "CHAIN_ID", "RES_ID", "X", "Y", "Z", "OCCUPANCY", 
                    "BETAFACTOR","CHARGE", "ELEMENT"]
    columsNums = [(0, 6), (6, 11), (11, 17), (17, 21), (21, 22), (22, 26), 
                  (26, 38), (38, 46), (46, 54), (54, 60), (60, 70), (70, 77), (77, 79)]

    # read pdbqt file        
    data = []
    with open(pdbqtFile, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                record = [line[start:end].strip() for start, end in columsNums]
                data.append(record)
    df = pd.DataFrame(data,columns=pdbqtColumns)
    # set appropriate types for elements in dataframe
    df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
    df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
    return df

#########################################################################################################################

def pdb_to_pdbqt(name, pdbFile, outDir, util24Dir, mglToolsDir,jobType,flexRes=None):
    env = os.environ.copy()
    env["PYTHONPATH"] = mglToolsDir
    os.chdir(outDir)
    prepReceptorPy = p.join(util24Dir, "prepare_receptor4.py")
    prepligandPy = p.join(util24Dir,"prepare_ligand4.py")
    prepFlexReceptorPy = p.join(util24Dir,"prepare_flexreceptor4.py")

    if jobType == "rigid":
        protPdbqt = p.join(outDir,"{}.pdbqt".format(name))
        call(["python2.7",prepReceptorPy,"-r",pdbFile,"-o",protPdbqt],
             env=env,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return protPdbqt

    elif jobType == "flex":
        if flexRes == None:
            #print(f"--X-->\tNo Flexible residues supplied for {name}")
            return
        rigidPdbqt = p.join(outDir,f"{name}_rigid.pdbqt")
        flexPdbqt = p.join(outDir,f"{name}_flex.pdbqt")
        call(["python2.7",prepFlexReceptorPy,"-r",pdbFile,"-s",flexRes,"-g",rigidPdbqt,"-x",flexPdbqt],
             env=env,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return rigidPdbqt, flexPdbqt
    
    elif jobType == "ligand":
        ligandPdbqt = p.join(outDir,f"{name}.pdbqt")
        call(["python2.7",prepligandPy,"-l",pdbFile,"-o",ligandPdbqt],
             env=env,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return ligandPdbqt