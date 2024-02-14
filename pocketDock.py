#########################################################################################################################
## import basic libraries
import os
import sys
import os.path as p
import pandas as pd
import multiprocessing as mp
import argpass
from tqdm import tqdm
## pocketDock modules
from pdbUtils import *
from modules_pocketDock import *
#########################################################################################################################
# get inputs
def read_inputs():
    # create an argpass parser, read config file, snip off ".py" if on the end of file
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()
    configName=args.config
    configName = p.splitext(configName)[0]

    # add config to PYTHONPATH
    cwd = os.getcwd()
    configPath = p.join(cwd,configName)
    sys.path.append(configPath)
    # import config file and run input function to return variables
    try:
        config_module = __import__(configName)
        (protDir, ligandDir, outDir, mglToolsDir, util24Dir,
             ligandOrdersCsv, modelSelectionMode, maxFlexRes,nCoresPerRun) = config_module.inputs()
        return (protDir, ligandDir, outDir, mglToolsDir, util24Dir,
             ligandOrdersCsv, modelSelectionMode, maxFlexRes,nCoresPerRun)
    except ImportError:
        print(f"Error: Can't to import module '{configName}'. Make sure the input exists!")
        print("HOPE IS THE FIRST STEP ON THE ROAD TO DISAPPOINTMENT")
        exit()

#########################################################################################################################
#########################################################################################################################
def main():
    (protDir, ligandDir, outDir, mglToolsDir, util24Dir,
             ligandOrdersCsv, modelSelectionMode, maxFlexRes,nCoresPerRun) = read_inputs()
    # disable copy warnings
    pd.set_option('mode.chained_assignment', None)

    exhausiveness, numModes = choose_model_selection_mode(modelSelectionMode)
    if not exhausiveness:
        print("Options for modelSelectionMode are \"best\", \"broad\", or \"balenced\"")
        exit()

    # make outDir
    os.makedirs(outDir,exist_ok=True)    
    # read ligand orders csv file into a dictionary
    ordersDf = pd.read_csv(ligandOrdersCsv)
    ordersDf["ID"] = ordersDf["ID"].astype(str)
    ordersDict = ordersDf.set_index('ID')['Ligand'].to_dict()



    pdbFiles =[]
    # loop through receptor PDB files
    for fileName in os.listdir(protDir):
        # Extract file name, skip if not a PDB file
        fileData = p.splitext(fileName)
        if not fileData[1] == ".pdb":
            continue
        pdbFiles.append(fileName)

    # run_paralell(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,
    #              mglToolsDir, exhausiveness, numModes,maxFlexRes,nCoresPerRun)

    run_serial(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,
     mglToolsDir, exhausiveness, numModes,maxFlexRes,nCoresPerRun)
#########################################################################################################################
#########################################################################################################################
def run_serial(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,
               mglToolsDir, exhausiveness, numModes, maxFlexRes,nCoresPerRun):
    # for testing 
    for fileName in pdbFiles:
        docking_protocol(fileName,protDir, ligandDir,outDir,ordersDict,util24Dir,
                         mglToolsDir,exhausiveness, numModes, maxFlexRes,nCoresPerRun)


#########################################################################################################################
def run_paralell(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,
                 mglToolsDir, exhausiveness, numModes,maxFlexRes,nCoresPerRun):
    num_cores = mp.cpu_count()
    paralellCores= round(num_cores / 2)

    with mp.Pool(processes=paralellCores) as pool:
        pool.starmap(docking_protocol,
                     tqdm([(fileName,protDir, ligandDir,outDir,ordersDict,util24Dir,
                            mglToolsDir,exhausiveness, numModes,maxFlexRes,nCoresPerRun) for fileName in pdbFiles],
                                total=len(pdbFiles)))


    # Close the pool to release resources
    pool.close()
    pool.join()


#########################################################################################################################

def docking_protocol(fileName,protDir, ligandDir,outDir,ordersDict,util24Dir,
                     mglToolsDir,exhausiveness, numModes, maxFlexRes,nCoresPerRun):
    # set up run directory and output key variables
    protName, protPdb, ligandPdb, ligandName, runDir = set_up_directory(fileName=fileName,
                                                                            protDir=protDir,
                                                                            ligandDir = ligandDir,
                                                                            outDir=outDir,
                                                                            ordersDict=ordersDict)  
    # Use fpocket to identify largest pocket, return center of pocket as [X,Y,Z] coords and Fpocket residues
    boxCenter, pocketResidues       =   run_fpocket(name=protName,
                                                        runDir=runDir,
                                                        pdbFile=protPdb)

    flexibeResidues =                   select_flexible_residues(protName=protName,
                                                                    protPdb=protPdb,
                                                                    flexResList=pocketResidues,
                                                                    maxFlexRes=maxFlexRes)

    protPdbqt                      =   pdb_to_pdbqt(name = protName,
                                                        pdbFile=protPdb,
                                                        outDir = runDir,
                                                        util24Dir = util24Dir,
                                                        mglToolsDir  = mglToolsDir,
                                                        jobType = "rigid")

    # Convert protein PDB to rigid and flexible PDBQT files
    rigidPdbqt, flexPdbqt                        =   pdb_to_pdbqt(name = protName,
                                                        pdbFile=protPdbqt,
                                                        outDir = runDir,
                                                        util24Dir = util24Dir,
                                                        mglToolsDir  = mglToolsDir,
                                                        jobType = "flex",
                                                        flexRes=flexibeResidues)

    # Convert ligand PDB to PDBQT files
    ligandPdbqt                     =   pdb_to_pdbqt(name = ligandName,
                                                        pdbFile=ligandPdb,
                                                        outDir = runDir,
                                                        util24Dir = util24Dir,
                                                        mglToolsDir  = mglToolsDir,
                                                        jobType = "ligand")

    # Write a config file for vina
    vinaConfig, dockedPdbqt         =   write_vina_config(outDir = runDir,
                                                            receptorPdbqt = rigidPdbqt,
                                                            flexPdbqt = flexPdbqt,
                                                            flex=True,
                                                            ligandPdbqt = ligandPdbqt,
                                                            boxCenter = boxCenter,
                                                            boxSize = 30,
                                                            exhaustiveness=exhausiveness,
                                                            numModes=numModes,
                                                            cpus=str(nCoresPerRun))

    # Run vina docking
    run_vina(outDir = runDir,
                configFile = vinaConfig)
    # split docking results PDBQT file into separate PDB files
    process_vina_results(outDir = runDir,
                            receptorPdbqt = protPdb,
                            dockedPdbqt = dockedPdbqt)


#########################################################################################################################

main()