import os
from os import path as p
import pdbUtils
from shutil import copy, rmtree
import pandas as pd

def clean_up_manager(cleanUpInfo, runDir, topDir):
    if cleanUpInfo["collateOutputs"]:
        collate_pdbs(runDir=runDir, outDir=topDir)
    if cleanUpInfo["genDockingReport"]:
        gen_docking_report(runDir=runDir, outDir=topDir)
    if cleanUpInfo["removeRunDirs"]:
        rmtree(runDir)
#########################################################################
def collate_pdbs(outDir, runDir):
    collatedDir = p.join(outDir,"collated_docked_pdbs")
    os.makedirs(collatedDir,exist_ok=True)
    pdbDir = p.join(runDir,"final_docked_pdbs")
    if not p.isdir(pdbDir):
        return
    for file in os.listdir(pdbDir):
        if not p.splitext(file)[1] == ".pdb":
            continue
        pdbFile = p.join(pdbDir,file)
        copy(pdbFile, p.join(collatedDir, file))
############################################################
def gen_docking_report(outDir, runDir):
    dirName = p.basename(runDir)
    index = 0
    energyDict = {}
    bindingPosePdbqt = p.join(runDir,"binding_poses.pdbqt")
    if not p.isfile(bindingPosePdbqt):
        return
    with open(bindingPosePdbqt,"r") as f:
        for line in f.readlines():
            if line.startswith("MODEL"):
                modelNumber = line.split()[1]
            if "VINA RESULT:" in line:
                bindingEnergy = float(line.split()[3])
                tmpDict = {"ID" : dirName,
                        "Binding Mode": modelNumber,
                        "Binding Energy": bindingEnergy}
                energyDict.update({index:tmpDict})
                index += 1
    reportDf = pd.DataFrame(energyDict).T
    print(reportDf)
    reportCsv = p.join(outDir,f"docking_report_{dirName}.csv")
    reportDf.to_csv(reportCsv, mode="w", index=False)
############################################################
def merge_report_dfs(outDir):
    dfsToConcat = []
    for file in os.listdir(outDir):
        if not file.startswith("docking_report"):
            continue
        csvFile = p.join(outDir,file)
        df = pd.read_csv(csvFile)
        dfsToConcat.append(df)
        os.remove(csvFile)
    reportDf = pd.concat(dfsToConcat,axis = 0, ignore_index=True)
    reportCsv = p.join(outDir, "docking_energies.csv")
    reportDf.to_csv(reportCsv,index=False)    