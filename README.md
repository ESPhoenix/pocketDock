# pocketDock
Modified Vina docking procedure, with pocket detection - fully automated
# PocketDock Installation Guide

## Step 1: Create Python Environment

Create a Python environment using Conda:

```bash
conda create -n pocketDock37 python=3.7
```

## Step 2: Install MGLTools

Download MGLTools from [here](https://ccsb.scripps.edu/mgltools/downloads/).

```bash
cd /home/{username}/bin
tar -xvzf mgltools_Linux-x86_64_1.5.7.tar.gz
cd mgltools_Linux-x86_64_1.5.7
chmod +x install.sh
./install.sh
```

## Step 3: Install Vina

Download AutoDock Vina from [here](https://vina.scripps.edu/downloads/).

```bash
cd /home/{username}/bin/
tar -xzf autodock_vina_1_1_2_linux_x86.tgz
cd autodock_vina_1_1_2_linux_x8/bin
chmod +x vina
```

Add Python2 to your PATH in .bashrc (or .bash_profile):

```bash
export PATH="$PATH:/usr/bin/python2.7"
source ~/.bashrc
```

Install Numpy for Python 2:

```bash
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
python2.7 get-pip.py
python2.7 -m pip install numpy
```

## Step 4: Install Python3 Libraries

Install required Python3 libraries using pip:

```bash
pip3 install pytest-shutil
pip3 install pandas
conda config --add channels conda-forge
conda install fpocket
pip3 install argpass
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
