import os
from os import path as p
import pdbUtils
from shutil import copy, rmtree
import pandas as pd
import multiprocessing as mp
import glob

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


########################################################################################
def merge_results(outDir, nCpus):
    csvFiles = glob.glob(os.path.join(outDir, "docking_report_*.csv"))
    csvBatches = batch_list(csvFiles,nCpus)
    outCsvs = [p.join(outDir,f"docking_energies_batch_{i}.csv") for i in range(1,len(csvBatches)+1)]

    with mp.Pool(processes=nCpus) as pool:
        jobs = []
        for csvFiles, outCsv in zip(csvBatches, outCsvs):
            job = pool.apply_async(merge_csv_files,(csvFiles,outCsv))
            jobs.append(job)
        for job in jobs:
            job.get()
    finalMergedCsv = p.join(outDir, "multiCaveFeatures.csv")
    merge_csv_files(outCsvs, finalMergedCsv)

########################################################################################
def merge_csv_files(csvFiles, outFile):
    chunk_size = min((10, len(csvFiles)))  # Adjust this based on available memory
    for i in range(0, len(csvFiles), chunk_size):
        chunkFiles = csvFiles[i:i + chunk_size]
        dfs = [pd.read_csv(csvFile) for csvFile in chunkFiles]
        merged_df = pd.concat(dfs, axis=0)
        merged_df.to_csv(outFile, mode='a', header=not os.path.exists(outFile), index=True)
    remainingFiles = csvFiles[i+chunk_size:]
    if not len(remainingFiles) == 0:
        dfs = [pd.read_csv(csvFile) for csvFile in remainingFiles]
        merged_df = pd.concat(dfs, axis=0) 
        merged_df.to_csv(outFile, mode='a', header=not os.path.exists(outFile), index=True)

    for csvFile in csvFiles:
        os.remove(csvFile)
###############################################################################################
def batch_list(inputList, nBatches):
    inLength = len(inputList)
    batchSize = inLength // nBatches
    remainder = inLength % nBatches

    batches = [inputList[i * batchSize + min(i, remainder):(i + 1) * batchSize + min(i + 1, remainder)]
                for i in range(nBatches)]
    return batches
