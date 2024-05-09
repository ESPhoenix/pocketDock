#########################################################################################################################
## import basic libraries
import os
import sys
import os.path as p
import pandas as pd
import multiprocessing as mp
import argpass
from tqdm import tqdm
import yaml
import time
## pocketDock modules
from pdbUtils import *
from modules_pocketDock import *
from cleanup_pocketDock import clean_up_manager, merge_report_dfs
#########################################################################################################################
# get inputs
def read_inputs():
    ## create an argpass parser, read config file, 
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    configFile=args.config
    ## Read config.yaml into a dictionary
    with open(configFile,"r") as yamlFile:
        config = yaml.safe_load(yamlFile) 

    return config
#########################################################################################################################
#########################################################################################################################
def main():
    print("POCKET DOCK TIME!")
    # read config from yaml file
    config = read_inputs()
    # extract needed paths
    pathInfo = config["pathInfo"]
    ligandDir = pathInfo["ligandDir"]
    outDir = pathInfo["outDir"]
    dockingOrders = pathInfo["dockingOrdersYaml"]
    with open(dockingOrders,"r") as yamlFile:
        dockingOrders = yaml.safe_load(yamlFile) 

    # disable copy warnings
    pd.set_option('mode.chained_assignment', None)
    # make outDir
    os.makedirs(outDir,exist_ok=True)    
    # pre-prepare ligand pdbqt files
    print("making ligand pdbqt files")
    gen_ligand_pdbqts(dockingOrders, ligandDir)
    # gererate a seqence of docking runs 
    print("Made ligand pdbqt files")

    if config["dockingInfo"]["totalCpus"] == 1:
        run_serial(config, dockingOrders)
    else:
        run_paralell(config, dockingOrders)

    # collect output files into single directory
    if config["cleanUpInfo"]["removeRunDirs"]:
        merge_report_dfs(outDir=outDir)
#########################################################################################################################
#########################################################################################################################
def run_serial(config, dockingOrders):#
    print("Running Serial")
    # for testing 
    for dockingOrder in dockingOrders:
        print(dockingOrder)
        docking_protocol(config, dockingOrder)


#########################################################################################################################
def run_paralell(config, dockingOrders):    
    print("running Parallel")
    totalCpus = config["dockingInfo"]["totalCpus"]
    cpusPerRun = config["dockingInfo"]["nCoresPerRun"]
    paralellCores = totalCpus // cpusPerRun

    with mp.Pool(processes=paralellCores) as pool:
        pool.starmap(docking_protocol,
                     tqdm([(config, dockingOrder) for dockingOrder in dockingOrders],
                                total=len(dockingOrders)))

#########################################################################################################################
def docking_protocol(config, dockingOrder):
    ## unpack config
    pathInfo = config["pathInfo"]
    outDir      = pathInfo["outDir"]

    ## unpack  docking parameters
    dockingInfo = config["dockingInfo"]
    maxFlexRes      = dockingInfo["maxFlexRes"]

    ## unpack general info
    generalInfo = config["generalInfo"]
    fpocketInputs = generalInfo["fpocketInputs"]

    pocketTag = None
    if fpocketInputs:
        pocketTag = dockingOrder["pocketTag"]


    # set up run directory and output key variables
    protName, protPdb, ligPdbqts, runDir = set_up_directory(outDir=outDir,
                                                            pathInfo=pathInfo,
                                                            dockingOrder=dockingOrder)  

    if  fpocketInputs:
        # get boxCenter and pocketResidues from pocketTag in docking instructions and fpocket pdb 
        pocketTag = dockingOrder["pocketTag"]
        # else:
        boxCenter, pocketResidues       = get_box_from_fpocket_inputs(pdbFile = protPdb,
                                                                    pocketTag = pocketTag)
        # remove fpocket pockets from pdb file
        protPdb = remove_fpockets(protPdb = protPdb,
                                   outDir = runDir,
                                    protName=protName)
  
    else:
        boxCenter, pocketResidues       =   run_fpocket(runDir=runDir,
                                                            pdbFile=protPdb)
    if maxFlexRes > 0:
        flexibeResidues =                   select_flexible_residues(protName=protName,
                                                                        protPdb=protPdb,
                                                                        flexResList=pocketResidues,
                                                                        maxFlexRes=maxFlexRes)
        
        rigidPdbqt, flexPdbqt           = gen_flex_pdbqts(protPdb = protPdb,
                                                flexibeResidues = flexibeResidues,
                                                outDir = runDir)

        # Write a config file for vina
        vinaConfig, dockedPdbqt         =   write_vina_config(outDir = runDir,
                                                                receptorPdbqt = rigidPdbqt,
                                                                flexPdbqt = flexPdbqt,
                                                                flex=True,
                                                                boxCenter = boxCenter,
                                                                boxSize = 30,
                                                                dockingInfo=dockingInfo)
    else:
        rigidPdbqt = pdb_to_pdbqt(inPdb=protPdb,
                                  outDir=runDir,
                                  jobType="rigid")
        

            # Write a config file for vina
        vinaConfig, dockedPdbqt         =   write_vina_config(outDir = runDir,
                                                                receptorPdbqt = rigidPdbqt,
                                                                flex=False,
                                                                boxCenter = boxCenter,
                                                                boxSize = 30,
                                                                dockingInfo=dockingInfo)
                                                            
    # Run vina docking
    run_vina(outDir = runDir,
                configFile = vinaConfig,
                ligPdbqts = ligPdbqts)
    # split docking results PDBQT file into separate PDB files
    dockedPdbs = process_vina_results(outDir = runDir,
                            receptorPdbqt = protPdb,
                            dockedPdbqt = dockedPdbqt,
                            dockingOrder = dockingOrder)
    cleanUpInfo = config["cleanUpInfo"]
    clean_up_manager(cleanUpInfo=cleanUpInfo,
                     runDir=runDir, 
                     topDir=outDir)

#########################################################################################################################
startTime = time.time()
main()
endTime = time.time()
elapsedTime = endTime - startTime
with open('time_output.txt', 'w') as f:  # Write elapsed time to file
    f.write(f"The script took {elapsedTime} seconds to run.")