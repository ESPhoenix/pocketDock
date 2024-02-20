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

Create a Python script for the configuration file, e.g., `config_file.py`. Fill in the required paths:

```python
def inputs():
    protDir = "path/to/pdb/files"
    ligandDir = "path/to/ligand/files"
    outDir = "path/to/desired/output/files"
    mglToolsDir = "/home/{username}/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs"
    util24Dir = "/home/{username}/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"
    vinaExe = "/home/{username}/bin/vina/vina_1.2.5_linux_x86_64"
    ligandOrdersCsv = "path/to/ligands_for_binding.csv"

    return protDir, ligandDir, outDir, mglToolsDir, util24Dir, vinaExe, ligandOrdersCsv
```

## Step 6: Run PocketDock

Run the PocketDock script with the provided configuration file:

```bash
python pocketDock.py --config config_file.py
```

**Note:** Replace `{username}` and update file paths accordingly in the above commands. Make sure to adjust permissions and paths based on your system configuration.
