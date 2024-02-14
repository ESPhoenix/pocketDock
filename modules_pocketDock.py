#########################################################################################################################
def choose_model_selection_mode(modelSelectionMode):
    if modelSelectionMode == "best":
        exhaustiveness = 16
        numModes = 1
    elif modelSelectionMode == "broad":
        exhaustiveness = 16
        numModes = 10
    elif modelSelectionMode == "balenced":
        exhaustiveness = 8
        numModes = 8
    else:
        return None, None
    return exhaustiveness, numModes
#########################################################################################################################
def process_vina_results(outDir,dockedPdbqt,receptorPdbqt,flex=False):
    # read output pdbqt file into a list of dataframes
    dockingDfList = read_docking_results(dockedPdbqt)
    receptorDf = pdbqt2df(receptorPdbqt)
    
    splice_docking_results(dockingDfList, receptorDf,outDir)
#########################################################################################################################
def splice_docking_results(dockingDfList, receptorDf, outDir):
    finalPdbDir = p.join(outDir,"final_docked_pdbs")
    os.makedirs(finalPdbDir,exist_ok=True)
    ## loop over each pose in dockingDfList
    for poseNumber, dockedDf in zip(range(1,len(dockingDfList)+1),dockingDfList):
        ## find  max chain ID in receptor, set ligand to one more than that
        lastChainIdInProt = receptorDf.iloc[-1]["CHAIN_ID"]
        ligandChainId = chr((ord(lastChainIdInProt) - ord('A') + 1) % 26 + ord('A'))
        dockedDf.loc[dockedDf["RES_ID"] == 0,"CHAIN_ID"] = ligandChainId
        ## find max residue number, set ligand residue number to one more than that
        lastResidueIdInProt = receptorDf.iloc[-1]["RES_ID"]
        ligandResidueId = int(lastResidueIdInProt)+1
        dockedDf.loc[dockedDf["RES_ID"] == 0,"RES_ID"] = ligandResidueId
        ## chage HETATM to ATOM for dockedDf
        dockedDf.loc[:,"ATOM"] = "ATOM"

        # Concat docked and rigid DFs togeter - this is in a weird order
        wholeDisorderedDf = pd.concat([dockedDf,receptorDf],axis=0)
        # get a list of unique residue Ids
        uniqueResidues = sorted(pd.unique(wholeDisorderedDf["RES_ID"]).tolist())
        # get one df per residue
        orderedResidues = []
        for residueNum in uniqueResidues:
            residueDf = wholeDisorderedDf[wholeDisorderedDf["RES_ID"]==residueNum]
            orderedResidues.append(residueDf)
        # concat into correct order
        wholeDf = pd.concat(orderedResidues)
        # re-do atom numbers
        wholeDf.loc[:,"ATOM_ID"] = range(1,len(wholeDf)+1)
        # save as pdb file
        saveFile = p.join(finalPdbDir, f"docked_pose_{str(poseNumber)}.pdb")
        df2Pdb(df=wholeDf,outFile=saveFile)
#########################################################################################################################
def read_docking_results(dockedPdbqt):
    # remove ROOT/BRANCH
    pdbqtColumns    =   ["ATOM","ATOM_ID", "ATOM_NAME", "RES_NAME",
                    "CHAIN_ID", "RES_ID", "X", "Y", "Z", "OCCUPANCY", 
                    "BETAFACTOR","CHARGE", "ELEMENT"]
    columsNums = [(0, 6), (6, 11), (11, 17), (17, 21), (21, 22), (22, 26), 
                  (26, 38), (38, 46), (46, 54), (54, 60), (60, 70), (70, 77), (77, 79)]
    # read pdbqt file into multiple dataframes
    dockingDfList =[]
    data = []
    # read output PDBQT file into set of dataframes
    with open(dockedPdbqt, 'r') as file:            
        for line in file:   
            if line.startswith("MODEL"):        # Each binding pose starts with "MODEL"
                if data == []:                  # skip 1st "MODEL"
                    continue
                df = pd.DataFrame(data,columns=pdbqtColumns)
                df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
                df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
                dockingDfList.append(df)
                data = []
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                record = [line[start:end].strip() for start, end in columsNums]
                data.append(record)
    # deal with last entry in pdbqtfile
    df = pd.DataFrame(data,columns=pdbqtColumns)
    df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
    df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
    dockingDfList.append(df)

    return dockingDfList
#########################################################################################################################
def run_vina(outDir,configFile):
    logFile = p.join(outDir,"vina_docking.log")
    with open(logFile,"a") as logFile:
        call(["vina","--config",configFile],stdout=logFile)
#########################################################################################################################
## writes a config file for a Vina docking simulation
def write_vina_config(outDir,receptorPdbqt,ligandPdbqt,boxCenter,boxSize,flexPdbqt=None,
                        exhaustiveness = 16, numModes = 10, cpus=2, energyRange = 5, seed = 42, flex=False):
    vinaConfigFile=p.join(outDir,f"vina_conf.txt")
    with open(vinaConfigFile,"w") as outFile:
        if flex:
            outFile.write(f"receptor = {receptorPdbqt}\n")
            outFile.write(f"flex = {flexPdbqt}\n\n")
            outFile.write(f"ligand = {ligandPdbqt}\n\n")
        else:
            outFile.write(f"receptor = {receptorPdbqt}\n")
            outFile.write(f"ligand = {ligandPdbqt}\n\n")

        outFile.write(f"center_x = {str(boxCenter[0])}\n")
        outFile.write(f"center_y = {str(boxCenter[1])}\n")
        outFile.write(f"center_z = {str(boxCenter[2])}\n\n")

        outFile.write(f"size_x = {str(boxSize)}\n")
        outFile.write(f"size_y = {str(boxSize)}\n")
        outFile.write(f"size_z = {str(boxSize)}\n\n")

        outFile.write(f"exhaustiveness = {str(exhaustiveness)}\n")
        outFile.write(f"num_modes = {str(numModes)}\n")
        outFile.write(f"energy_range = {str(energyRange)}\n\n")
        outFile.write(f"seed = {str(seed)}\n\n")

        if not flex:
            dockedPdbqt  =   p.join(outDir,f"binding_poses.pdbqt")
            outFile.write(f"out = {dockedPdbqt}\n")
        else:
            dockedPdbqt  =   p.join(outDir,f"binding_poses.pdbqt")
            outFile.write(f"out = {dockedPdbqt}\n")
        outFile.write(f"cpu = {cpus}")

        return vinaConfigFile, dockedPdbqt
#########################################################################################################################
def select_flexible_residues(protName,protPdb,flexResList,maxFlexRes):

    priorityFlexRes = ["TRP", "TYR", "PHE",         ## AROMATICS
                       "ARG","LYS","HIS",           ## LARGE, POSITIVE
                       "GLN","GLU","ASN","ASP",     ## NEGATIVE OR H-BONDING
                       "MET","LEU","ILE","VAL",     ## LARGE, HYDROPHOBIC
                       "SER","CYS","THR",           ## SMALL, PROTIC
                       "ALA","PRO","GLY"]           ## SMALL, NO ROTATION IN SIDE CHAIN - NEVER GOING TO USE THESE
    # load protein PDB as a DF
    protDf = pdb2df(protPdb)
    # reduce DF down to unique residues in flexResList (created from binding pocket)
    flexDf = protDf[protDf["RES_ID"].isin(flexResList)]
    flexDf.drop_duplicates(subset=["RES_ID"], inplace=True)
    # sort by priorityFlexRes, get up to the maximum number of flexible residues 
    flexDf.loc[:,'RES_NAME'] = pd.Categorical(flexDf['RES_NAME'], categories=priorityFlexRes, ordered=True)
    flexDf = flexDf.sort_values(by="RES_NAME")
    flexDf = flexDf.head(maxFlexRes)
        
    # get a list of unique CHAIN IDs
    uniqueChains = flexDf.drop_duplicates(subset=["CHAIN_ID"])["CHAIN_ID"].tolist()
    # generate MGLTools compatable flexible residues
    flexResidues=[]
    for chain in uniqueChains:
        chainBlock=[]
        for _, row in flexDf.iterrows():
            if row["CHAIN_ID"]==chain:
                resName = row["RES_NAME"]
                resSeq = row["RES_ID"]
                selection=f"{resName}{resSeq}"
                chainBlock.append(selection)
        chainBlock="_".join(chainBlock)
        chainBlock=f"{protName}:{chain}:{chainBlock}"
        flexResidues.append(chainBlock)
    
    flexResidues=",".join(flexResidues)
    return flexResidues    
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
    
    #########################################################################################################################
def set_up_directory(fileName,protDir,ligandDir,outDir,ordersDict):
    fileData = p.splitext(fileName)
    protPdb = p.join(protDir,fileName)
    protName = fileData[0]
    runDir = p.join(outDir,protName)
    os.makedirs(runDir,exist_ok=True)
    copy(protPdb,runDir)
    protPdb = p.join(runDir,f"{protName}.pdb")

    # identify ligand 
    ligandName = ordersDict[protName]
    ligandPdb = p.join(ligandDir,f"{ligandName}.pdb")
    copy(ligandPdb,runDir)
    ligandPdb = p.join(runDir,f"{ligandName}.pdb")

    return protName, protPdb, ligandPdb, ligandName, runDir


#########################################################################################################################

def run_fpocket(name,runDir,pdbFile):
   # print("----->\tRunning Fpocket!")
    os.chdir(runDir)
    minSphereSize = "3.0"
    maxSphereSize = "6.0"
    call(["fpocket","-f",pdbFile,"-m",minSphereSize,"-M",maxSphereSize],
         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    fpocketOutDir = p.join(runDir,f"{name}_out","pockets")
    ## ASSUMPTION == LARGEST POCKET IS OUR BINDING POCKET ## Not really true!
    largestPocketPdb = p.join(fpocketOutDir,"pocket1_atm.pdb")
    ## ERROR Handling
    if not p.isfile(largestPocketPdb):
        #print("--X-->\tFpocket Failed!")
        return
    largestPocketDf = pdb2df(largestPocketPdb)

    boxCenter = [largestPocketDf["X"].mean(), largestPocketDf["Y"].mean(),largestPocketDf["Z"].mean()]
    pocketResidues = largestPocketDf["RES_ID"].unique().tolist()

    #print("----->\tFpocket Success!")
    return boxCenter, pocketResidues