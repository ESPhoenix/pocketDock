def inputs():
    # protDir = location of your protein PDB files
    protDir = "/home/esp/dataset_generation/flavin_dataset/04_3d_structure_generation/Alphafold/01_50_percent_similarity/fixed_folds/All_folds"
    # ligandDir = location of your ligand PDB files
    ligandDir = "/home/esp/dataset_generation/flavin_dataset/06_cofactor_complexes/ligands"
    # where do you want output files to be put (will make this dir if it doesn't already exist)
    outDir = "/home/esp/dataset_generation/flavin_dataset/06_cofactor_complexes/outputs"
    # location of your MGLToolsPckgs directory
    mglToolsDir = "/home/esp/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs"
    # location of AutoDockTools/Utilities24 firectory from MGLTools
    util24Dir = "/home/esp/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"
    # location of Vina executable
    vinaExe = "/home/esp/bin/vina/vina_1.2.5_linux_x86_64"
    # location of CSV with columns "ID" and "Cofactor" for each enzyme-substrate complex you want to make
    ligandOrdersCsv = "/home/esp/dataset_generation/flavin_dataset/05_cofactor_ID/cofactor_predictions.csv"
    # Tunes exhaustiveness and num_modes parameters, options: "best", "broad"."balenced"
    modelSelectionMode = "best"     

    return protDir, ligandDir, outDir, mglToolsDir, util24Dir, vinaExe, ligandOrdersCsv, modelSelectionMode