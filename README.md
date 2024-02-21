# pocketDock
Modified Vina docking procedure, with pocket detection - fully automated
# PocketDock Installation Guide

## Step 1: Create Python Environment

Create a Python environment using Conda:

```bash
conda create -n pocketDock37 python=3.7
```


## Step 2: Install Python3 Libraries

Install required Python3 libraries using pip and conda:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install openbabel
conda install fpocket
conda install vina
pip install pyyaml
pip install argpass
pip install tqdm
pip install pytest-shutil
pip install pandas
pip install pyarrow
```

## Step 5: Create Config File

Create a Yaml file for the configuration file, e.g., `config_file.yaml`. Fill in the required paths, parameters for docking and options for post-docking clean-up:

```yaml
pathInfo:
  protDir: "/home/esp/scriptDevelopment/pocketDock/receptors"
  ligandDir: "/home/esp/scriptDevelopment/pocketDock/ligands"
  outDir: "/home/esp/scriptDevelopment/pocketDock/outputs"
  ligandOrdersCsv: "/home/esp/scriptDevelopment/pocketDock/cofactor_predictions.csv"
dockingInfo:
  maxFlexRes: 5
  nCoresPerRun: 2 
  totalCpus: 24
  exhaustiveness: 16
  numModes: 1
cleanUpInfo:
  collateOutputs: True
  genDockingReport: True
  removeRunDirs: True
```

## Step 6: Run voidDock

Run the PocketDock script with the provided configuration file:

```bash
python voidDock.py --config config_file.py
```

