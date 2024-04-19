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
## pocketDock modules
from pdbUtils import *
from modules_pocketDock import *
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
    # read config from yaml file
    config = read_inputs()
    # extract needed paths
    pathInfo = config["pathInfo"]
    protDir = pathInfo["protDir"]
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
    gen_ligand_pdbqts(dockingOrders, ligandDir)
    # gererate a seqence of docking runs 
    dockingSequence = gen_docking_sequence(dockingOrders, protDir, ligandDir)

    if config["dockingInfo"]["totalCpus"] == 1:
        run_serial(config, dockingSequence)
    else:
        run_paralell(config, dockingSequence)

    # collect output files into single directory
    clean_up(config["cleanUpInfo"], outDir)
#########################################################################################################################
#########################################################################################################################
def run_serial(config, dockingSequence):
    # for testing 
    for dockId in dockingSequence:
        docking_protocol(config, dockDetails = dockingSequence[dockId])


#########################################################################################################################
def run_paralell(config, dockingSequence):

    paralellCores= config["dockingInfo"]["totalCpus"]

    dockingDetailsList = []
    for dockId in dockingSequence:
        dockingDetailsList.append(dockingSequence[dockId])
    with mp.Pool(processes=paralellCores) as pool:
        pool.starmap(docking_protocol,
                     tqdm([(config, dockDetails) for dockDetails in dockingDetailsList],
                                total=len(dockingSequence)))

#########################################################################################################################
def docking_protocol(config, dockDetails):
    print(dockDetails)

    ## unpack config
    pathInfo = config["pathInfo"]
    outDir      = pathInfo["outDir"]

    ## unpack  docking parameters
    dockingInfo = config["dockingInfo"]
    maxFlexRes      = dockingInfo["maxFlexRes"]

    ## unpack general info
    generalInfo = config["generalInfo"]
    fpocketInputs = generalInfo["fpocketInputs"]

    # ## deal with multiple ligands being docked in order
    # if len(dockDetails["ligPdb"]) > 1:
    #     ligPdb = dockDetails["ligPdb"][0]
    #     remainingLigPdbs = dockDetails["ligPdb"][1:]
    # else:
    #     ligPdb = dockDetails["ligPdb"][0]
    #     remainingLigPdbs = False
    # dockDetails["ligPdb"] = ligPdb


    # set up run directory and output key variables
    protName, protPdb, ligPdbqts, runDir = set_up_directory(outDir=outDir,
                                                            pathInfo=pathInfo,
                                                            dockDetails=dockDetails)  

    # Use fpocket to identify largest pocket, return center of pocket as [X,Y,Z] coords and Fpocket residues
    if  fpocketInputs:
        # get boxCenter and pocketResidues from pocketTag in docking instructions and fpocket pdb 
        pocketTag = dockDetails["pocketTag"]
        # if "boxCenter" in dockDetails:
        #     boxCenter = dockDetails["boxCenter"]
        #     pocketResidues = dockDetails["pocketResidues"]
        # else:
        boxCenter, pocketResidues       = get_box_from_fpocket_inputs(pdbFile = protPdb,
                                                                    pocketTag = pocketTag)
        # remove fpocket pockets from pdb file
        protPdb = remove_fpockets(protPdb = protPdb,
                                   outDir = runDir,
                                    protName=protName)
    else:
        if "boxCenter" in dockDetails:
            boxCenter = dockDetails["boxCenter"]
            pocketResidues = dockDetails["pocketResidues"]
        else:
            boxCenter, pocketResidues       =   run_fpocket(runDir=runDir,
                                                                pdbFile=protPdb)

    # if "flexibleResidues" in dockDetails:
    #     flexibeResidues = dockDetails["flexibleResidues"]
    # else:
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
                                                            
    # Run vina docking
    run_vina(outDir = runDir,
                configFile = vinaConfig,
                ligPdbqts = ligPdbqts)
    # split docking results PDBQT file into separate PDB files
    dockedPdbs = process_vina_results(outDir = runDir,
                            receptorPdbqt = protPdb,
                            dockedPdbqt = dockedPdbqt,
                            dockDetails = dockDetails)
    
    # ## use outputs for more ligands to be docked:
    # if remainingLigPdbs:
    #     newRunId = dockDetails["runId"] + 1
    #     newDockingDetails = {"protPdb":dockedPdbs[0],
    #                         "ligPdb":remainingLigPdbs,
    #                         "pocketResidues": pocketResidues,
    #                         "flexibleResidues": flexibeResidues,
    #                         "boxCenter": boxCenter,
    #                         "runId": newRunId,
    #                         "pocketTag": pocketTag}
    #     docking_protocol(config=config, dockDetails=newDockingDetails)


#########################################################################################################################

main()