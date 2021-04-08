from scipy.stats import kendalltau
import pprint,os,math
import numpy as np
from keras.constraints import max_norm
from numpy.random import seed
from keras.models import Sequential
from keras.models import load_model
from keras.layers import Dense, Activation,Dropout
from keras.optimizers import SGD,Adam
from keras.wrappers.scikit_learn import KerasRegressor
from keras.callbacks import ReduceLROnPlateau
import os
from scipy import stats
from argparse import RawTextHelpFormatter
from keras.layers import merge, Conv1D, MaxPooling2D, Input, Dense, Flatten, Reshape, RepeatVector, LSTM, embeddings, MaxPooling1D, LocallyConnected1D, Lambda
import keras.layers
from keras import regularizers
from keras.models import Model, Sequential, load_model
from keras.callbacks import ModelCheckpoint  
from keras.utils.np_utils import to_categorical
from keras.layers import Activation, Dropout
from keras.optimizers import SGD, Adam,Nadam
from keras.metrics import top_k_categorical_accuracy
from keras import backend as K
def readexpt(exptfile):
    def readchain(chainlist):
        Sumpro_dict={}
        for chain in chainlist:
            chainAB=chain[5:];pro_dict={}
            for i in range(len(chainAB)):
                if chainAB[i]=='_':
                    #ProA=chainAB[:i];ProB=chainAB[i+1:]
                    pro_dict['ProA']=chainAB[:i]
                    pro_dict['ProB']=chainAB[i+1:]
            Sumpro_dict[chain[:4]]=pro_dict
        return Sumpro_dict
    PDBsum=[];PDBID=[];exptdict=dict();lines=open(exptfile,'r').readlines()
    for line in lines:
        line=str(line).strip().split()
        PDBID.append(str(line[0])[:4])
        PDBsum.append(str(line[0]))
    Pdbid=np.unique(PDBID)
    PDBsum=np.unique(PDBsum)
    Pro_dict=readchain(PDBsum)
    pprint.pprint(Pro_dict)
    for pdbname in Pdbid:
        son_dict=dict()
        for line in lines:
            line=str(line).strip().split()
            if str(pdbname)==str(line[0][:4]):
                #son_dict=dict()
                son_dict[str(line[1])]=str(line[3])
        exptdict[str(pdbname)]=son_dict
    #pprint.pprint(exptdict)
    return lines,exptdict,Pro_dict
def findfiles(path,label):
    files=[];filenames=os.listdir(path)
    for name in filenames:
        if os.path.splitext(name)[1]==str(label):
            files.append(name)
    return files
def creatinput_matrixc(file,pdblist):
    def findvdw_r(restype):
        res_atom={'ALA':['CT'],'A':['CT'],'ARG':['C8','C8','C8','C8','CA','N2','N2','N2'],'R':['C8','C8','C8','C8','CA','N2','N2','N2'],'ASP':['2C','O2','O2','CO'],'D':['2C','O2','O2','CO'],'CYS':['2C','SH'],'C':['2C','SH'],'GLN':['2C','2C','N','C','O'],'Q':['2C','2C','N','C','O'],'GLU':['2C','2C','CO','O2','O2'],'E':['2C','2C','CO','O2','O2'],'LYS':['C8','C8','C8','N3'],'K':['C8','C8','C8','N3'],'HID':['CT','CC','CV','CR','NB'],'HIE':['CT','CC','CV','CR','NB','NA'],'HIS':['CT','CC','CV','CR','NB','NA'],'HIP':['CT','CC','CV','CR','NA','NA'],'H':['CT','CC','CV','CR','NB','NA'],'ILE':['3C','2C','CT','CT'],'I':['3C','2C','CT','CT'],'LEU':['2C','3C','CT','CT'],'L':['2C','3C','CT','CT'],'MET':['2C','2C','CT','S'],'M':['2C','2C','CT','S'],'PHE':['CT','CA','CA','CA','CA','CA','CA'],'F':['CT','CA','CA','CA','CA','CA','CA'],'PRO':['CT','CT','CT'],'P':['CT','CT','CT'],'SER':['2C','OH'],'S':['2C','OH'],'THR':['3C','CT','OH'],'T':['3C','CT','OH'],'TRP':['CT','CA','CA','CA','CA','CB','CN','CW','C','NA'],'W':['CT','CA','CA','CA','CA','CB','CN','CW','C','NA'],'TYR':['CT','CA','CA','CA','CA','CA','C','OH'],'Y':['CT','CA','CA','CA','CA','CA','C','OH'],'VAL':['3C','CT','CT'],'V':['3C','CT','CT']}
        lines=open('vdw_r.txt','r').readlines()
        vdw_r={}
        for line in lines:
            line=str(line).strip().split()
            vdw_r[str(line[0])]=str(line[2])
        vdw_l=np.zeros(10)
        if restype in res_atom.keys():
            atomlist=res_atom[restype]
            for i in range(len(atomlist)):
                vdw_l[i]=float(vdw_r[atomlist[i]])
        return vdw_l
    lines=open(file,'r').readlines()
    res_dict={'H':[0,0,0,0,0],'HIS':[0,0,0,0,0],'HIE':[0,0,0,0,0],'HID':[0,0,0,0,0],'HIP':[0,0,0,0,0],'R':[0,0,0,0,1],'ARG':[0,0,0,0,1],'D':[0,0,0,1,0],'ASP':[0,0,0,1,0],'E':[0,0,0,1,1],'GLU':[0,0,0,1,1],\
              'LYS':[0,0,1,0,0],'K':[0,0,1,0,0],'C':[0,0,1,0,1],'CYS':[0,0,1,0,1],'Q':[0,0,1,1,0],'GLN':[0,0,1,1,0],'N':[0,0,1,1,1],'ASN':[0,0,1,1,1],'S':[0,1,0,0,0],'SER':[0,1,0,0,0],'T':[0,1,0,0,1],'THR':[0,1,0,0,1],\
              'Y':[0,1,0,1,0],'TYR':[0,1,0,1,0],'L':[0,1,0,1,1],'LEU':[0,1,0,1,1],'M':[0,1,1,0,0],'MET':[0,1,1,0,0],'F':[0,1,1,0,1],'PHE':[0,1,1,0,1],'I':[0,1,1,1,0],'ILE':[0,1,1,1,0],\
              'W':[0,1,1,1,1],'TRP':[0,1,1,1,1],'V':[1,0,0,0,0],'VAL':[1,0,0,0,0],'G':[1,0,0,0,1],'GLY':[1,0,0,0,1],'A':[1,0,0,1,0],'ALA':[1,0,0,1,0],'P':[1,0,0,1,1],'PRO':[1,0,0,1,1]}
    seconds_dict={'None':[0,0,0,0],'C':[0,0,0,1],'E':[0,0,1,0],'T':[0,0,1,1],'H':[0,1,0,0],'G':[0,1,0,1],'B':[0,1,1,0],'b':[0,1,1,1],'I':[1,0,0,0]}
    aminotype={'I':[0,0],'ILE':[0,0],'F':[0,0],'PHE':[0,0],'L':[0,0],'LEU':[0,0],'W':[0,0],'TRP':[0,0],'A':[0,0],'ALA':[0,0],'M':[0,0],'MET':[0,0],'P':[0,0],'PRO':[0,0],'V':[0,0],'VAL':[0,0],\
               'C':[0,1],'CYS':[0,1],'N':[0,1],'ASN':[0,1],'G':[0,1],'GLY':[0,1],'S':[0,1],'SER':[0,1],'Q':[0,1],'GLN':[0,1],'Y':[0,1],'TYR':[0,1],'T':[0,1],'THR':[0,1],\
               'H':[1,0],'HIS':[1,0],'HIP':[1,0],'HIE':[1,0],'HID':[1,0],'R':[1,0],'ARG':[1,0],'K':[1,0],'LYS':[1,0],'D':[1,1],'ASP':[1,1],'E':[1,1],'GLU':[1,1]}
    Cluster=[];mutadict={};mutationnumber=0
    for i in range(len(lines)):
        line=str(lines[i])
        line=line.strip().split()       
        if line[0]=='PDBID:' and line[1][:4]+'.pdb' in pdblist:
            #f1.write(str(line) +'\n')
            mutationnumber += 1
            cluster=[];row=[];count=1
            expt=float(line[7]);row.append(expt)
            pre_res=res_dict[line[3]];row+=pre_res
            pre_type=aminotype[line[3]];row+=pre_type
            now_res1=res_dict[line[5]];row+=now_res1
            now_type=aminotype[line[5]];row+=now_type
            second1=seconds_dict[line[12]];row+=second1
            sinphi1=math.sin((float(line[14])/180)*np.pi);row.append(sinphi1)
            cosphi1=math.cos((float(line[14])/180)*np.pi);row.append(cosphi1)
            sinpsi1=math.sin((float(line[16])/180)*np.pi);row.append(sinpsi1)
            cospsi1=math.cos((float(line[16])/180)*np.pi);row.append(cospsi1)
            sasa1=float(line[18]);row.append(sasa1/50.0)
            sasa2=float(line[20])-sasa1;row.append(sasa2/50.0)
            pssmwt=line[24];row.append(float(pssmwt)/10.0)
            pssmdiff=line[26];row.append(float(pssmdiff)/10.0)
            while True:
                line1=str(lines[i+count]).strip().split()
                if line1[0]=='Around_res:':
                    #Around_res:  A3 CAdis:  3.789 Structure:  C Phi:  -138.94 Psi:  148.75 SASA:  88.800 SASApart : 88.800 Hbonds:  0 protein: 1
                    row_oth=[];sumrow=[]
                    llCa=str(lines[i+count+2]).strip().split()
                    res_around=str(llCa[3])[-3:]
                    row_oth+=res_dict[res_around]
                    res_type=aminotype[str(llCa[3])[-3:]];row_oth+=res_type
                    Xca=float(llCa[6]);Yca=float(llCa[7]);Zca=float(llCa[8])
                    llN=str(lines[i+count+1]).strip().split()
                    Xn=float(llN[6]);Yn=float(llN[7]);Zn=float(llN[8])
                    llC=str(lines[i+count+3]).strip().split()
                    Xc=float(llC[6]);Yc=float(llC[7]);Zc=float(llC[8])

                    Xn=Xn-Xca;Yn=Yn-Yca;Zn=Zn-Zca
                    Xc=Xc-Xca;Yc=Yc-Yca;Zc=Zc-Zca
                    lenca_n=np.sqrt(Xn**2+Yn**2+Zn**2)
                    lenca_c=np.sqrt(Xc**2+Yc**2+Zc**2)
                    lenca_ca=np.sqrt(Xca**2+Yca**2+Zca**2)
                    Xn=Xn/lenca_n;Yn=Yn/lenca_n;Zn=Zn/lenca_n####vector from Ca to N
                    Xc=Xc/lenca_c;Yc=Yc/lenca_c;Zc=Zc/lenca_c####vector from Ca to C
                    Xca=Xca/lenca_ca;Yca=Yca/lenca_ca;Zca=Zca/lenca_ca####vector from Ca to Ca of the center residue

                    CaCadis=float(line1[3]);row_oth.append(float(CaCadis)/10.0)
                    second2=seconds_dict[line1[5]];row_oth+=second2
                    sinphi2=math.sin((float(line1[7])/180)*np.pi);row_oth.append(sinphi2)
                    cosphi2=math.cos((float(line1[7])/180)*np.pi);row_oth.append(cosphi2)
                    sinpsi2=math.sin((float(line1[9])/180)*np.pi);row_oth.append(sinpsi2)
                    cospsi2=math.cos((float(line1[9])/180)*np.pi);row_oth.append(cospsi2)
                    sasa1=float(line1[11]);row_oth.append(sasa1/50.0)
                    sasa2=float(line1[14])-sasa1;row_oth.append(sasa2/50.0)
                    row_oth.append(Xn);row_oth.append(Yn);row_oth.append(Zn)
                    row_oth.append(Xc);row_oth.append(Yc);row_oth.append(Zc)
                    row_oth.append(Xca);row_oth.append(Yca);row_oth.append(Zca)
                    count=count+1
                    sumrow=row+row_oth
                    sumrow=list(map(float,sumrow))
                    cluster.append(sumrow)
                elif line1[0]=='PDBID:' or int(i+count)>=len(lines)-1:
                    break
                else:
                    count += 1
            cluster=np.array(cluster)
          
            mutadict[mutationnumber]=cluster
    ClusterNo=25
    for i in range(ClusterNo):
        inputx=[]
        if i==0:
            RealY=[]
            for j in range(mutationnumber):
                RealY.append(mutadict[j+1][0][0])
        for j in range(mutationnumber):
            row=[]
            row=mutadict[j+1][i][1:]
            inputx.append(row)
        inputx=np.array(inputx)
        Cluster.append(inputx)
    return Cluster,ClusterNo,RealY

def pred(pdbname):
    train_t,dictumber,trainY_t=creatinput_matrixc('Newallaroundres_test.txt',[pdbname])
    model = load_model('ddg_ppi.h5')
    ppred=model.predict(train_t,verbose=1)
    print('predict:',ppred)
    return ppred

#pred('1a22.pdb')
