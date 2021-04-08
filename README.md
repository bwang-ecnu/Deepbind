# Deepbind
Machine learning model for predicting the effection of point mutation to protein-protein binding free energy 

Requirements

tensorflow > 1.4 or 2.x

keras 2.4.3

naccess

python 3

Usage
1. Install naccess (http://www.bioinf.manchester.ac.uk/naccess)

tar -zxvf naccess.tar.gz

Then change the second line in 'naccess' file: set set EXE_PATH = /home/bwang/github/Naccess(your path of Naccess)

2. Make predictoin
python predict.py 1a22.pdb(change with your padfile)
