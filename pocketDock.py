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
    pathInfo = config["pathInfo"]
    dockingInfo = config["dockingInfo"]
    return pathInfo,dockingInfo
#########################################################################################################################
#########################################################################################################################
def main():
    pathInfo, dockingInfo  = read_inputs()
    protDir = pathInfo["protDir"]
    ligandDir = pathInfo["ligandDir"]
    outDir = pathInfo["outDir"]
    ligandOrdersCsv = pathInfo["ligandOrdersCsv"]
    # disable copy warnings
    pd.set_option('mode.chained_assignment', None)
    # make outDir
    os.makedirs(outDir,exist_ok=True)    

    # gererate a seqence of docking runs 
    dockingSequence = gen_docking_sequence(ligandOrdersCsv, protDir, ligandDir)
    run_paralell(pathInfo, dockingInfo, dockingSequence)

    #run_serial(pathInfo, dockingInfo, dockingSequence)
#########################################################################################################################
#########################################################################################################################
def run_serial(pathInfo, dockingInfo, dockingSequence):
    # for testing 
    for dockId in dockingSequence:
        docking_protocol(pathInfo, dockingInfo, dockDetails = dockingSequence[dockId])


#########################################################################################################################
def run_paralell(pathInfo, dockingInfo, dockingSequence):
    paralellCores= dockingInfo["totalCpus"]

    dockingDetailsList = []
    for dockId in dockingSequence:
        dockingDetailsList.append(dockingSequence[dockId])
    with mp.Pool(processes=paralellCores) as pool:
        pool.starmap(docking_protocol,
                     tqdm([(pathInfo, dockingInfo, dockDetails) for dockDetails in dockingDetailsList],
                                total=len(dockingSequence)))

#########################################################################################################################

def docking_protocol(pathInfo, dockingInfo, dockDetails):
    ## unpack path Info
    outDir      = pathInfo["outDir"]
    mglToolsDir = pathInfo["mglToolsDir"]
    util24Dir   = pathInfo["util24Dir"]

    ## unpack general docking parameters
    maxFlexRes      = dockingInfo["maxFlexRes"]
    exhaustiveness  = dockingInfo["exhaustiveness"]
    numModes        = dockingInfo["numModes"]
    nCoresPerRun    = dockingInfo["nCoresPerRun"]
    
    # set up run directory and output key variables
    protName, protPdb, ligPdb, ligandName, runDir = set_up_directory(outDir=outDir,
                                                                            dockDetails=dockDetails)  
    # Use fpocket to identify largest pocket, return center of pocket as [X,Y,Z] coords and Fpocket residues
    boxCenter, pocketResidues       =   run_fpocket(runDir=runDir,
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
    ligPdbqt                     =   pdb_to_pdbqt(name = ligandName,
                                                        pdbFile=ligPdb,
                                                        outDir = runDir,
                                                        util24Dir = util24Dir,
                                                        mglToolsDir  = mglToolsDir,
                                                        jobType = "ligand")

    # Write a config file for vina
    vinaConfig, dockedPdbqt         =   write_vina_config(outDir = runDir,
                                                            receptorPdbqt = rigidPdbqt,
                                                            flexPdbqt = flexPdbqt,
                                                            flex=True,
                                                            ligPdbqt = ligPdbqt,
                                                            boxCenter = boxCenter,
                                                            boxSize = 30,
                                                            exhaustiveness=exhaustiveness,
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