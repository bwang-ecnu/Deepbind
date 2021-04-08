import os,re,math
import numpy as np

#Collect the expt results
def readexpt(exptfile):
    def readchain(chainlist):
        Sumpro_dict={}
        for chain in chainlist:
            chainAB=chain[5:];pro_dict={}
            for i in range(len(chainAB)):
                if chainAB[i]=='_':
                    pro_dict['ProA']=chainAB[:i]
                    pro_dict['ProB']=chainAB[i+1:]
            Sumpro_dict[chain[:4]]=pro_dict
        return Sumpro_dict
    PDBsum=[];PDBID=[];exptdict=dict();exptdict2=dict();lines=open(exptfile,'r').readlines()
    for line in lines:
        line=str(line).strip().split()
        PDBID.append(str(line[0])[:4])
        PDBsum.append(str(line[0]))
    Pdbid=np.unique(PDBID)
    PDBsum=np.unique(PDBsum)
    Pro_dict=readchain(PDBsum)
    for pdbname in Pdbid:
        son_dict=dict()
        son_dict2=dict()
        for line in lines:
            line=str(line).strip().split()
            if str(pdbname)==str(line[0][:4]):
                son_dict[str(line[1])]=str(line[3])
                if str(line[1])[:1]!='C':
                    son_dict2[str(line[1])]=str(line[3])
        exptdict[str(pdbname)]=son_dict
        exptdict2[str(pdbname)]=son_dict2
    print(exptdict)
    print(Pro_dict)
    return lines,exptdict,exptdict2,Pro_dict
#Find all PDB files
def findpdbfile(path,label):
    files=[];filenames=os.listdir(path)
    for name in filenames:
        if os.path.splitext(name)[1]==str(label):
            files.append(name)
    return files
#Check to make no reduant PDB file
def clean(path,Lexpt,Lpdbfile):
    Pdbid=[]
    os.system('mkdir Newpdbs')
    path1='./Newpdbs/'
    for line in Lexpt:
        line=line.strip().split()
        Pdbid.append(line[0])
    PDBID=np.unique(Pdbid)
    aminolist=['HIS','HID','HIP','ARG','LYS','ILE','PHE','LEU','TRP','ALA','MET','PRO','CYX','CYS','ASN','VAL'\
               ,'GLY','SER','GLN','TYR','THR','GLU','ASP']
    for pid in PDBID:
        for name in Lpdbfile:
            if str(pid)[:4]==str(name)[:4]:
                Ca=[];N=[];C=[];O=[]
                Com=[];Main_atom=['CA','N','C','O']
                lines1=open(str(path)+'/'+name,'r').readlines()
                f=open(path1+str(name),'w')
                for i in range(len(lines1)-2):
                    l1=str(lines1[i]).strip().split();l2=str(lines1[i+1]).strip().split();l3=str(lines1[i+2]).strip().split()
                    ll1=str(lines1[i]);ll2=str(lines1[i+1])
                    if ('ATOM' in lines1[i]):
                        if (('1.00' in lines1[i])):pass
                        elif (str(l1[2]).strip() in Main_atom) and (str(l2[2]).strip()) and (str(l1[2])==str(l2[2])) and (float(ll1[56:61])==0.50):
                            Com.append(lines1[i+1])
                        elif (str(l1[2]) in Main_atom) and  (float(ll1[56:61])==0.25) and (float(ll2[56:61])==0.25) and (str(l1[2])==str(l2[2])):
                            testline=str(lines1[i+2])
                            if float(testline[56:61])==0.25:
                                Com.append(lines1[i]);Com.append(lines1[i+1]);Com.append(lines1[i+2])
                            elif float(testline[56:61])==0.50:
                                Com.append(lines1[i]);Com.append(lines1[i+1])
                        elif (str(l1[2]) in Main_atom) and (float(ll1[56:61])==0.33) and (str(l1[2])==str(l2[2])==str(l3[2])):
                            Com.append(lines1[i+2]);Com.append(lines1[i+1])
                        elif (str(l1[2]) in Main_atom) and (str(l1[2])==str(l2[2])) and (float(ll1[56:61])>float(ll2[56:61])):
                            Com.append(lines1[i+1])
                        elif (str(l1[2]) in Main_atom) and (str(l1[2])==str(l2[2])) and (float(ll1[56:61])<float(ll2[56:61])):
                            Com.append(lines1[i])
                for line in lines1:
                    lls=str(line).strip().split()
                    if 'ATOM' in lls:
                        if (line not in Com) and (lls[3][-3:] in aminolist) :
                            f.write(str(line))
                f.close()
                lines1=open(str(path)+'/'+name,'r').readlines()
                for line in lines1: 
                    if 'ATOM' in line and line not in Com:
                        l=line.strip().split()
                        if 'CA'==str(l[2]).strip():Ca.append(line)
                        elif 'N'==str(l[2]).strip():N.append(line)
                        elif 'C'==str(l[2]).strip():C.append(line)
                        elif 'O'==str(l[2]).strip():O.append(line)
                if len(Ca)!=len(N)!=len(C)!=len(O):
                    print(pid);print(len(Ca));print(len(N))
                elif len(Ca)==len(N)==len(C)==len(O):print('Clean '+str(name)+ ' Finish!')
####function to get mutation pdbfile
def rosetta(Exptdict):
    os.system('mkdir Rosetta_pdb')
    for pdbname in Exptdict.keys():
        pdbfile=pdbname+'.pdb'
        for residno in Exptdict[pdbname].keys():  
            no=str(residno[2:-1])
            preres=str(residno[:1])
            mutares=str(residno[-1])
            chain=str(residno[1:2])
            if preres!='C':
                f=open('resfile.txt','w')
                f.write('NATRO' +'\n')
                f.write('start' +'\n')
                f.write('    '+str(no)+' '+str(chain)+' PIKAA '+str(mutares))
                f.close()
                Rosetta_Muta.ChangeA_R(pdbfile,1)
                Rosetta_Muta.Rosetta_PDB('AtoR.pdb')
                Rosetta_Muta.ChangeA_R('AtoR_0001.pdb',0)
                os.system('mv  RtoA.pdb  Rosetta_pdb/'+str(pdbname)+'_'+str(chain)+str(no)+str(mutares)+'.pdb')
                os.system('rm AtoR.pdb');os.system('rm AtoR_0001.pdb')

#sort function to get the around residues and their coord
def sorttest(Calist,chainID,NoID,Xa,Ya,Za):
    sort_dis=dict();sort_xyz=dict()
    for line in Calist:
        line=str(line)
        if (str(line[21:22])+str(line[22:26]).strip())!=(str(chainID)+str(NoID)):
            aroundNo=str(line[21:22])+str(line[22:26]).strip();X1=float(line[30:38]);Y1=float(line[38:46]);Z1=float(line[46:54])
            distance=np.sqrt((X1-Xa)**2+(Y1-Ya)**2+(Z1-Za)**2)
            sort_dis[aroundNo]=distance;sort_xyz[distance]=line[30:38]+' '+line[38:46]+' '+line[46:54]
    sort_dis=sorted(sort_dis.items(),key=lambda x: x[1])
    sort_xyz=sorted(sort_xyz.items(),key=lambda x: x[0])
    return sort_dis,sort_xyz
####function to find secondary structure
def findsecondstructure(choice,file,CHain,NO):
    if choice==1:
        f=open(file,'r')
        lines=f.readlines()
    elif choice==2:
        f=open('Muta_stride/'+str(file),'r')
        lines=f.readlines()
    STRUC='None';PHI=0;PSI=0
    for line in lines:
        line=str(line).strip().split()
        if 'ASG'==line[0] and str(CHain)==line[2] and str(NO)==line[3]:
            STRUC=line[5];PHI=line[7];PSI=line[8]        
    f.close()
    return STRUC,float(PHI),float(PSI)
####function to find SASA
def findsasa(choice,file,CHain,NO):
    if choice==1:
        f=open(file,'r');lines=f.readlines();sasa=0
        for line in lines:
            line=str(line).strip().split()
            if 'RES'==line[0] and str(CHain)==line[2] and str(NO)==line[3]:
                sasa=line[5]
        f.close()
    elif choice==2:
        f=open('Mutasasa/'+str(file),'r');lines=f.readlines();sasa=0
        for line in lines:
            line=str(line).strip().split()
            if 'RES'==line[0] and str(CHain)==line[2] and str(NO)==line[3]:
                sasa=line[5]
        f.close()
    print('CHAIN,NO',CHain,NO)
    print('sasa:',sasa)
    return float(sasa)
####function to find the number of H bonds
def findHbonds(choice,file,Chain,No,Chain_oth,No_oth):
    if choice==1:
        f=open(file,'r')
        lines=f.readlines()
    elif choice==2:
        f=open('MutaHbplus/'+str(file),'r')
        lines=f.readlines()
    Hbonds=0
    def changelen(Number):
        if len(Number)==1:
            Number='000'+Number
        elif len(Number)==2:
            Number='00'+Number
        elif len(Number)==3:
            Number='0'+Number
        return Number
    No=changelen(str(No));No_oth=changelen(str(No_oth))
    Inmutares=Chain+No;Incontres=Chain_oth+No_oth    
    for line in lines:
        line=str(line).strip().split()
        mutares=line[0][:5];contres=line[2][:5]
        if (Inmutares==mutares and Incontres==contres) or (Inmutares==contres and Incontres==mutares):
            Hbonds += 1
    f.close()
    return Hbonds
####function to find the PSSM_wt and PSSM_diff
def findpssm(pdbname,chain,no,preres,nowres):
    def findfile(path,label):
        files=[];filenames=os.listdir(path)
        for name in filenames:
            if os.path.splitext(name)[1]==str(label):
                files.append(name)
        return files
    pssmdict={'A':2,'R':3,'N':4,'D':5,'C':6,'Q':7,'E':8,'G':9,'H':10,'I':11,'L':12,'K':13,'M':14,'F':15,'P':16,'S':17,'T':18,'W':19,'Y':20,'V':21}
    pssmfiles=findfile('./pssmresults','.pssm')
    prescore=0;scorediff=0
    for pssmfile in pssmfiles:
        if str(pdbname)==pssmfile[:4] and str(chain) in str(pssmfile[5:-5]):
            lines=open('pssmresults/'+str(pssmfile),'r').readlines()
            lines=lines[3:]
            for line in lines:
                if str(line).strip()!='':
                    line=str(line).strip().split()
                    if str(line[0])==str(no):
                        prescore=line[pssmdict[preres]]
                        nowscore=line[pssmdict[nowres]]
                        scorediff=float(nowscore)-float(prescore)
    return prescore,scorediff
def findpssm_oth(pdbname,chain,no,res):
    def findfile(path,label):
        files=[];filenames=os.listdir(path)
        for name in filenames:
            if os.path.splitext(name)[1]==str(label):
                files.append(name)
        return files
    pssmdict={'A':2,'R':3,'N':4,'D':5,'C':6,'Q':7,'E':8,'G':9,'H':10,'I':11,'L':12,'K':13,'M':14,'F':15,'P':16,'S':17,'T':18,'W':19,'Y':20,'V':21}
    pssmfiles=findfile('./pssmresults','.pssm')
    score=0
    for pssmfile in pssmfiles:
        if str(pdbname)==pssmfile[:4] and str(chain) in str(pssmfile[5:-5]):
            lines=open('pssmresults/'+str(pssmfile),'r').readlines()
            lines=lines[3:]
            for line in lines:
                if str(line).strip()!='':
                    line=str(line).strip().split()
                    if str(line[0])==str(no):
                        score=line[pssmdict[res]]
    return score
#find all PDB mute residue's and write the around residues' normalized coords to Allaroundres.txt
def find_localresidues(path,PDBfiles,PDBfiles2,Exptdict,path2,Exptdict2,Prochaindict):
    f=open('Allaroundres.txt','w')    
    for pdbid in PDBfiles:
        lines=open(path+'/'+pdbid,'r').readlines()
        CAlines=[line for line in lines if 'CA' in line and 'ATOM' in line]
        for pdbname in Exptdict.keys():
            if str(pdbid)[:4]==str(pdbname):
                muteno=[];pre_res=[];now_res=[];exp=[]
                print('PDBID:',pdbid)
                for residno in Exptdict[pdbname].keys():  
                    muteno.append(residno[1:-1])
                    pre_res.append(residno[:1])
                    now_res.append(residno[-1])
                    exp.append(Exptdict[pdbname][residno])
                for i in range(len(muteno)):
                    chain=str(muteno[i][:1]);No=str(muteno[i][1:])
                    preres=pre_res[i];nowres=now_res[i];expvalue=float(exp[i])
                    #print('Chain and No:',chain,No)
                    #pre_score,diff_score=findpssm(pdbname,chain,No,preres,nowres)
                    pre_score,diff_score=0,0
                    if chain in Prochaindict[pdbname]['ProA']:
                        protein=1
                        SASApart=findsasa(1,pdbname+'1.rsa',chain,No)
                    elif chain in Prochaindict[pdbname]['ProB']:
                        protein=2
                        SASApart=findsasa(1,pdbname+'2.rsa',chain,No)
                    Struc,PHi,PSi=findsecondstructure(1,pdbname+'.stride',chain,No)
                    SASA=findsasa(1,pdbname+'.rsa',chain,No)
                    f.write('{} {} {} {} {} {} {} {:.3f} {} {} {} {} {} {} {} {} {} {} {:.3f} {} {:.3f} {} {} {} {} {} {}'.format('PDBID: ',pdbid,'pre_res: ',\
                            preres,'now_res: ',nowres,'expt: ',expvalue,'ChainNo: ',chain,No,'Structure: ',\
                            Struc,'Phi: ',PHi,'Psi: ',PSi,'SASA: ',SASA,'SASApart: ',SASApart,'protein:',protein,'PSMMwt:',pre_score,'PSSMdiff:',diff_score) +'\n')
                    for line in CAlines:
                        line=str(line)
                        if line[21:22]==chain and str(line[22:26]).strip()==No:
                            xx=float(line[30:38].strip());yy=float(line[38:46].strip());zz=float(line[46:54].strip())
                    for line in lines:
                        line=str(line).strip()
                        if line[21:22]==chain and line[22:26].strip()==No:
                            xxx=float(line[30:38].strip());yyy=float(line[38:46].strip());zzz=float(line[46:54].strip())
                            f.write('{} {} {:.3f} {:.3f} {:.3f} {}'.format(line[:22],line[22:26],xxx,yyy,zzz,line[56:80]) +'\n')
                    Sort_resids,Sort_coor=sorttest(CAlines,chain,No,xx,yy,zz)
                    for i in range(len(Sort_resids[:25])):
                        chain_oth=str(Sort_resids[i][0])[:1];No_oth=str(Sort_resids[i][0])[1:]                        
                        if chain_oth in Prochaindict[pdbname]['ProA']:
                            protein=1
                            SASApart=findsasa(1,pdbname+'1.rsa',chain_oth,No_oth)
                        elif chain_oth in Prochaindict[pdbname]['ProB']:
                            protein=2
                            SASApart=findsasa(1,pdbname+'2.rsa',chain_oth,No_oth)                                                        
                        Struc,PHi,PSi=findsecondstructure(1,pdbname+'.stride',chain_oth,No_oth)
                        SASA=findsasa(1,pdbname+'.rsa',chain_oth,No_oth)
                        Hbonds=findHbonds(1,pdbname+'.hb2',chain,No,chain_oth,No_oth)
                        f.write('{} {} {} {:.3f} {} {} {} {} {} {} {} {:.3f} {} {:.3f} {} {} {} {}'.format('Around_res: ',Sort_resids[i][0],'CAdis: ',\
                                Sort_resids[i][1],'Structure: ',Struc,'Phi: ',PHi,'Psi: ',PSi,'SASA: ',SASA,'SASApart :',SASApart,'Hbonds: ',Hbonds,'protein:',protein)+'\n')
                        for line in lines:
                            line=str(line).strip()
                            if line[21:22]==chain_oth and str(line[22:26]).strip()==No_oth:
                                xxx=float(line[30:38].strip());yyy=float(line[38:46].strip());zzz=float(line[46:54].strip())
                                f.write('{} {} {:.3f} {:.3f} {:.3f} {}'.format(line[:22],line[22:26],xxx,yyy,zzz,line[56:80]) +'\n')                                          

  
    f.write('End')
    f.close()
####################
#translation matrix
def get_translationmatrix(vectora):
    #vectora=vectora/np.sqrt(np.sum(np.dot(vectora,vectora)))
    dx=float(vectora[0]);dy=float(vectora[1]);dz=float(vectora[2])
    trans_matrix=[]
    trans_matrix.append([1,0,0,0]);trans_matrix.append([0,1,0,0]);trans_matrix.append([0,0,1,0])
    trans_matrix.append([-dx,-dy,-dz,1])
    trans_matrix=np.array(trans_matrix)
    return trans_matrix
#rotation matrix
def get_rotatematrix(axis, origin, angle):
   angle = angle/180.0*3.1415926
   ccos = math.cos(angle)
   ssin = math.sin(angle)                                            
   a = np.array(axis)
   u, v, w = a / np.sqrt(np.sum(np.dot(a,a)))
   a, b, c = origin
   rot_matrix = []
   rot_matrix.append([u**2+(v**2+w**2)*ccos, u*v*(1-ccos)-w*ssin, u*w*(1-ccos)+v*ssin, \
           (a*(v**2+w**2)-u*(b*v+c*w))*(1-ccos)+(b*w-c*v)*ssin])
   rot_matrix.append([ u*v*(1-ccos)+w*ssin, v**2+(u**2+w**2)*ccos, v*w*(1-ccos)-u*ssin,\
       (b*(u**2+w**2)-v*(a*u+c*w))*(1-ccos)+(c*u-a*w)*ssin])
   rot_matrix.append([ u*w*(1-ccos)-v*ssin, v*w*(1-ccos)+u*ssin, w**2+(u**2+v**2)*ccos,\
           (c*(u**2+v**2)-w*(a*u+b*v))*(1-ccos)+(a*v-b*u)*ssin])
   rot_matrix.append([0, 0, 0, 1])
   rot_matrix = np.array(rot_matrix)
   return rot_matrix                                   
#orient the coord from Allaroundres.txt
def rotation(pre_coordfile):
    lines=open(pre_coordfile,'r').readlines()
    f=open('Newallaroundres_test.txt','w')
    #find the CA,N,C coord
    def getatomcoord(coordlist):
        Cacoor=[];Ncoor=[];Ccoor=[];Sumcoor=[]
        for line in coordlist:
            line=str(line).strip().split()
            if len(line)==12:
                Sumcoor.append(line[6].strip());Sumcoor.append(line[7].strip());Sumcoor.append(line[8].strip())
                if 'CA'==line[2].strip():
                    Cacoor.append(line[6].strip());Cacoor.append(line[7].strip());Cacoor.append(line[8].strip())
                elif 'N'==line[2].strip():
                    Ncoor.append(line[6].strip());Ncoor.append(line[7].strip());Ncoor.append(line[8].strip())
                elif 'C'==line[2].strip():
                    Ccoor.append(line[6].strip());Ccoor.append(line[7].strip());Ccoor.append(line[8].strip())
            elif len(line)==11:
                Sumcoor.append(line[5].strip());Sumcoor.append(line[6].strip());Sumcoor.append(line[7].strip())
                if 'CA'==line[2].strip():
                    Cacoor.append(line[5].strip());Cacoor.append(line[6].strip());Cacoor.append(line[7].strip())
                elif 'N'==line[2].strip():
                    Ncoor.append(line[5].strip());Ncoor.append(line[6].strip());Ncoor.append(line[7].strip())
                elif 'C'==line[2].strip():
                    Ccoor.append(line[5].strip());Ccoor.append(line[6].strip());Ccoor.append(line[7].strip())
        Cacoor=list(map(float,Cacoor));Ncoor=list(map(float,Ncoor));Ccoor=list(map(float,Ccoor))
        Cacoor=np.array(Cacoor,dtype=np.float32);Ncoor=np.array(Ncoor,dtype=np.float32);Ccoor=np.array(Ccoor,dtype=np.float32)
        Sumcoor=list(map(float,Sumcoor));Sumcoor=np.array(Sumcoor,dtype=np.float32);Sumcoor=Sumcoor.reshape(-1,3)
        return Cacoor,Ncoor,Ccoor,Sumcoor
    ####get the mute residue's coord and the translation and rotation matrixs
    for i in range(len(lines)):
        if 'PDBID:' in lines[i]:
            tempcoordlist=[];count=1
            f.write(str(lines[i]))
            while True:                
                if 'ATOM' not in lines[i+count]:
                    break
                else:
                    tempcoordlist.append(lines[i+count]);count+=1
            ####below to get the translation matrixs ,move Ca to x,y,z=(0,0,0)
            CAxyz,Nxyz_pre,Cxyz_pre,Sumxyz=getatomcoord(tempcoordlist)
            Mtrione=np.ones(Sumxyz.shape[0]);Sumxyz=np.c_[Sumxyz,Mtrione]            
            TransMatrix=get_translationmatrix(CAxyz)
            Maxt1=np.dot(Sumxyz,TransMatrix)
            ####below to get the rotation matrixs,align N to -x axis
            Nxyz=Maxt1[0,:3]
            Nxyz_len=np.sqrt(np.sum(np.dot(Nxyz,Nxyz)))
            Nxyz_norm=Nxyz/Nxyz_len
            xaxis=np.array([-1,0,0])
            cross=np.cross(Nxyz_norm,xaxis)
            cos=np.dot(Nxyz_norm,xaxis)
            angle_nx=np.arccos(cos)*180/np.pi
            cross1=np.cross(xaxis,cross)
            cos=np.dot(Nxyz_norm,cross1)
            angle1=np.arccos(cos)*180/np.pi
            if angle1<=90:
                rotangle=angle_nx
            else:
                rotangle=360-angle_nx
            RoatMatrix1=get_rotatematrix(cross,[0,0,0],rotangle)
            Maxt2=np.dot(Maxt1,RoatMatrix1.T)            
            ####below to align C to x-y plane and make CB at z>0 
            Cxyz=Maxt2[2,:3]
            cos=Cxyz[1]/np.sqrt(Cxyz[1]**2+Cxyz[2]**2)
            angle=np.arccos(cos)*180/np.pi
            if Cxyz[2]>0:
                angle=180-angle
            elif Cxyz[2]<=0:
                angle=180+angle
            RoatMatrix2=get_rotatematrix([1,0,0],[0,0,0],angle)           
            np.set_printoptions(precision=3,suppress=True)
            Maxt3=np.dot(Maxt2,RoatMatrix2.T)
            for j in range(int(Maxt3.shape[0])):
                line=str(tempcoordlist[j]).strip()
                f.write('{} {:.3f} {:.3f} {:.3f}'.format(line[:27]+'    ',Maxt3[j][0],Maxt3[j][1],Maxt3[j][2]) +'\n')               
            print('Sumxyz:',Sumxyz)
            print('Max3:',Maxt3)                               
            ####multiply other residues coord maxt with Translation and two Rotation Maxts
        elif 'Around_res:' in lines[i]:
            tempcoordlist=[];count=1
            f.write(str(lines[i]))
            while True:
                if 'ATOM' not in lines[i+count]:
                    break
                elif (i+count)<len(lines):
                    tempcoordlist.append(lines[i+count]);count+=1
            CAxyz,Nxyz_pre,Cxyz_pre,Sumxyz=getatomcoord(tempcoordlist)
            Matrione=np.ones(Sumxyz.shape[0]);Sumxyz=np.c_[Sumxyz,Matrione]
            Maxtros1=np.dot(Sumxyz,TransMatrix)
            Maxtros2=np.dot(Maxtros1,RoatMatrix1.T)
            Maxtros3=np.dot(Maxtros2,RoatMatrix2.T)
            print('Maxtros3:',Maxtros3)
            for j in range(int(Maxtros3.shape[0])):
                line=str(tempcoordlist[j]).strip()
                f.write('{} {:.3f} {:.3f} {:.3f}'.format(line[:27]+'    ',Maxtros3[j][0],Maxtros3[j][1],Maxtros3[j][2]) +'\n')
    f.write('End')
    f.close()
####use stride produce secondary structure files
def stride(choice,path,PdbList):
    if choice==1:
        for file in PdbList:
            os.system('./stride '+path+'/'+str(file)+' > '+str(file)[:-4]+'.stride')
    elif choice==2:
        os.system('rm -r Muta_stride')
        os.system('mkdir Muta_stride')
        for file in PdbList:
            os.system('./stride '+path+'/'+str(file)+' > '+'Muta_stride/'+str(file)[:-4]+'.stride')
####use Naccess produce SASA information files
def Naccess(chioce,path,PdbList,ProchainDict):
    ####sasa of protein AB
    for file in PdbList:
        os.system('./naccess '+path+'/'+str(file))
        ####divide the AB pdb file into A and B
        pdbname=str(file)[:4]
        writename=str(file)[:-4]
        f1=open(writename+'1.PDB','w')
        f2=open(writename+'2.PDB','w')
        proA=ProchainDict[pdbname]['ProA']
        proB=ProchainDict[pdbname]['ProB']
        lines=open(path+'/'+file,'r').readlines()
        for i in range(len(lines)):
            line=str(lines[i]).strip().split()
            if 'ATOM' in line:
                if str(line[4]) in str(proA):
                    #print('proA:',lines[i])
                    f1.write(str(lines[i]))
                elif str(line[4]) in str(proB):
                    #print('proB:',lines[i])
                    f2.write(str(lines[i]))            
        f1.close()
        f2.close()
        ####sasa of protein A
        os.system('./naccess '+writename+'1.PDB')
        ####sasa of protein B
        os.system('./naccess '+writename+'2.PDB')
        if chioce==2:
            os.system('mv '+pdbname+'_* Mutasasa')
####use Hbplus to calculate the number of h-bonds
def Hbplus(chioce,path,PdbList):
    for file in PdbList:
        os.system('./hbplus '+path+'/'+str(file))
        if chioce==2:
            pdbname=str(file)[:4]
            os.system('mv '+pdbname+'_* MutaHbplus')
####use python2.7 to fasten PDB files
def fasten(path,PdbList):
    for file in PdbList:
        pdbname=str(file)[:4]
        os.system('nohup python2.7 fasta.py '+path+'/'+str(file)+' '+pdbname+'.fasta  >  '+pdbname+'.log')
def fastenout(PdbList):
    sequence={}
    for file in PdbList:
        outfile=str(file)[:4]+'.log'
        lines=open(outfile,'r').readlines()
        for line in lines:
            if 'chaininfor' in line:
                line=line[12:50].strip()
                chains=''.join(re.findall(r'[A-Za-z]',line))
                if len(chains)>1:
                    sequence[str(file)[:4]]=chains
                    print(file)
                    print(chains) 
    print(sequence)
    
####main####                
if __name__=='__main__':
    print('welcome to collecting python program!')
    ##FIND THE INPUT DATA(EXPT AND PDB)##
    exptlist,exptDict,exptDict2,prochainDict=readexpt('no_reduan.txt')
    Pdblist=findpdbfile('./PDBs','.pdb')
    clean('./PDBs',exptlist,Pdblist)
    Newpdblist=findpdbfile('./Newpdbs','.pdb')
    Mutapdblist=[]
    print('collect the input information')
    ##COLLECT THE INPUT INFORMATION OF DATA SET AND WRITE TO A FILE##
    stride(1,'Newpdbs',Newpdblist)
    Naccess(1,'Newpdbs',Newpdblist,prochainDict)
    Hbplus(1,'Newpdbs',Newpdblist)
    find_localresidues('Newpdbs',Newpdblist,Mutapdblist,exptDict,'Rosetta_pdb',exptDict2,prochainDict)
    rotation('Allaroundres.txt')
