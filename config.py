def inputs():
    # protDir = location of your protein PDB files
    protDir = "/home/esp/scriptDevelopment/pocketDock/receptors"
    # ligandDir = location of your ligand PDB files
    ligandDir = "/home/esp/scriptDevelopment/pocketDock/ligands"
    # where do you want output files to be put (will make this dir if it doesn't already exist)
    outDir = "/home/esp/scriptDevelopment/pocketDock/outputs"
    # location of your MGLToolsPckgs directory
    mglToolsDir = "/home/esp/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs"
    # location of AutoDockTools/Utilities24 firectory from MGLTools
    util24Dir = "/home/esp/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"
    # location of CSV with columns "ID" and "Cofactor" for each enzyme-substrate complex you want to make
    ligandOrdersCsv = "/home/esp/scriptDevelopment/pocketDock/docking_commands.csv"
    # Tunes exhaustiveness and num_modes parameters, options: "best", "broad"."balenced"
    modelSelectionMode = "best"     
    # maximum number of flexible residues per docking simulation (time scales with this)
    maxFlexRes = 5
    # number of cores per docking simulation (use 1 for lots of receptors!)
    nCoresPerRun = 4

    return (protDir, ligandDir, outDir, mglToolsDir, util24Dir,
             ligandOrdersCsv, modelSelectionMode, maxFlexRes,nCoresPerRun)