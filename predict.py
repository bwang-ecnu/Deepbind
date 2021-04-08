import tensorflow.compat.v1 as tf
import keras
from keras.models import load_model
import numpy as np
import os
from ddg_pred import pred
import sys

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
model = load_model('ddg_ppi.h5')

def make_mutation(pdbfile,name):
    onelabel=['H','R','K','I','F','L','W','A','M','P','C','N','V','G','S','Q','Y','D','E','T']
    aminodict={'HIS':'H','HIP':'H','HIE':'H','ARG':'R','LYS':'K','ILE':'I','PHE':'F','LEU':'L','TRP':'W','ALA':'A','MET':'M','PRO':'P','CYS':'C','ASN':'N','VAL':'V','GLY':'G','SER':'S','GLN':'Q','TYR':'Y','ASP':'D','GLU':'E','THR':'T'}
    lines=open(pdbfile,'r').readlines()
    res_no=[];amino=[];chain=[]
    for line in lines:
        if 'ATOM' in line[:6]:
            line=line.strip().split()
            chain.append(line[4])
    chain=np.unique(np.array(chain))
    schain=''
    for i in chain:
        schain=schain+'_'+str(i)
    mutainfo=[]
    for line in lines:
        if 'ATOM' in line[:6]:
            line=line.strip().split()
            chain=line[4]
            resno=line[5];ami=aminodict[line[3]]
            for muta in onelabel:
                info=name[:-4]+schain+' '+ami+chain+resno+str(muta)+' none 0.0'
                if info not in res_no:
                    if str(muta)!=ami:
                        mutainfo.append(ami+chain+resno+str(muta))
                        res_no.append(info)    
    f=open('no_reduan.txt','w')
    for i in res_no:
        f.write(i+'\n')
    f.close()
    os.system('mkdir PDBs')
    os.system('cp '+pdbfile+' PDBs')   
    os.system('python prepare.py')    
    return mutainfo


def make(filename):
        mutainfo=make_mutation(filename,filename)
        result = pred(filename)
        predict=[]
        for i in range(len(result)):
            predict.append(round(float(result[i][0])*10,3))
        f=open('ddg_pred.txt','w')
        f.write('mutation ddg(kcal/mol) ddg=dgmutation-dgwildtype\n')
        for i in range(len(mutainfo)):
            f.write(str(mutainfo[i])+'   '+str(predict[i])+'\n') 
        f.close()

if __name__ == "__main__":
    pdb = sys.argv[1]
    make(pdb)
