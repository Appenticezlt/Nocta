#!/bin/bash
echo "Setting up environment for Nihilism  Project..."
conda create -n analysis python=3.10 -y
conda activate analysis
pip install -r requirements.txt

echo " Installing PyRosetta via pyrosetta-installer..."
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

echo "PyRosetta installation completed!"
echo “Environment has already sucessfully equipped！”