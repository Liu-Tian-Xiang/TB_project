import matplotlib.pyplot as plt
import numpy as np
import os

def getData(RunModes,bandbegin,bandend,LengthBegin,LengthEnd,dataDIR,Shift):
    data1=0
    aimDIR=dataDIR+str(bandbegin)+'_'+str(LengthBegin)
    aimFILE=dataDIR+str(bandbegin)+'_'+str(LengthBegin)+'/M'+str(RunModes)+'band'+str(0)+'.data'
    print (aimFILE)
    if os.path.isdir(aimDIR) and os.path.isfile(aimFILE):
        data1=np.genfromtxt(aimFILE)

    for j in range(1,int(Shift)):    
        aimFILE=dataDIR+str(bandbegin)+'_'+str(LengthBegin)+'/M'+str(RunModes)+'band'+str(j)+'.data'
        if os.path.isdir(aimDIR) and os.path.isfile(aimFILE):
            data2=np.genfromtxt(aimFILE)
            if data2.ndim!=1:
                data1=np.append(data1,data2,0)

    return data1

def getString(string,filename):
    with open(filename,'r') as f:
        for line in f:
            sline=line.split(' ')
            if string==sline[0]:
                n=sline[2]
                break
    f.close
    return n

def getList(string,filename):
    with open(filename,'r') as f:
        for line in f:
            sline=line.split(' ')
            if string==sline[0]:
                n=[int(x) for x in sline[2:]]
                break
    f.close
    return n

def getStr(string,filename):
    with open(filename,'r') as f:
        for line in f:
            sline=line.split(' ')
            if string==sline[0]:
                n=sline[2]
                break
    f.close
    return float(n)

def drawDiagrams(dataDIR,axes,DataFile,taglist,setcolor,titleTag):
    RunModes=getList('RunModes','./out.code')
    minDataTag=0
    maxDataTag=len(taglist)
    for pesudoNoc in range(minDataTag,maxDataTag):#number of curves
        noc=taglist[pesudoNoc]
        DividedFactor=getStr('DividedFactor',dataDIR+'0/out'+str(noc)+'.txt')
        totalblocks=int(2**(int(DividedFactor)))
        length=getStr('numprocs',dataDIR+'0/out'+str(noc)+'.txt')
        n=getStr('down',dataDIR+'0/out'+str(noc)+'.txt')
        m=getStr('up',dataDIR+'0/out'+str(noc)+'.txt')
        pnumber=int(getStr('pnumber',dataDIR+'0/out'+str(noc)+'.txt'))
        atom=getStr('AtomNum',dataDIR+'0/out'+str(noc)+'.txt')
        LengthBegin=getStr('LengthBegin',dataDIR+'0/out'+str(noc)+'.txt')
        LengthEnd=getStr('LengthEnd',dataDIR+'0/out'+str(noc)+'.txt')
        InAsLength=0
        if(LengthBegin==LengthEnd): InAsLength=LengthBegin
        bandbegin=getStr('bandbegin',dataDIR+'0/out'+str(noc)+'.txt')
        bandend=getStr('bandend',dataDIR+'0/out'+str(noc)+'.txt')
        AutoMatic=getStr('AutoMatic',dataDIR+'0/out'+str(noc)+'.txt')
        if(AutoMatic==0):
            bandbegin=0
            bandend=0
            LengthBegin=0
            LengthEnd=0
        DataTag=getStr('DataTag',dataDIR+'0/out'+str(noc)+'.txt')
        kpath=getString('kppath',dataDIR+'0/out'+str(noc)+'.txt')
        kvalue=getList('kvalue',dataDIR+'0/out'+str(noc)+'.txt')
        kpath=kpath[0:-1].split('->')
        kpos=kvalue[kpath.index('G')]
        kpos_in_id=0
        ExtendedRange=300
        NoteRange=1
    
        gappos=int(getStr('GapLevel'+str(int(DataTag)),dataDIR+'0/out'+str(noc)+'.txt'))+1
        rangv=gappos-1
        rangc=gappos
        kpath=getString('kppath',dataDIR+'0/out'+str(noc)+'.txt')
        kvalue=getList('kvalue',dataDIR+'0/out'+str(noc)+'.txt')
        kpath=kpath[0:-1].split('->')
        o=0
        for kv in kpath:
            if kv =='G':
                kpath[o]='$\Gamma$'
            o+=1
        for md in range(0,len(RunModes)):
            for w in range(int(LengthBegin),int(LengthEnd)+1):
                for l in range(int(bandbegin),int(bandend)+1):
                    GapLevel=gappos-1
                    print(GapLevel)
                    if GapLevel<ExtendedRange:
                        ExtendedRange=GapLevel
                    data1=getData(RunModes[md],l+int(DataTag),l+int(DataTag),w,w,dataDIR,length)
                    if setcolor!='none' and setcolor!='all':
                        line1,=plt.plot(data1[:,0],data1[:,rangv-ExtendedRange+1],label='Band gap['+str(data1[kpos,rangc]-data1[kpos,rangv])+'eV]',color=setcolor[pesudoNoc])
                    else:
                        line1,=plt.plot(data1[:,0],data1[:,rangv-ExtendedRange+1],label='Band gap['+str(data1[kpos,rangc]-data1[kpos,rangv])+'eV]')
                    clr=line1.get_color()
                    #for ii in range(1,2*ExtendedRange):
                    for ii in range(1,data1.shape[1]+ExtendedRange-rangv-1):
                        plt.plot(data1[:,0],data1[:,rangv-ExtendedRange+ii+1],color=clr)
                    area=(8)**2
                    vv=0
                    cc=0
                    for ii in range(0,NoteRange):
                        plt.scatter(data1[kpos,0],data1[kpos,rangv-NoteRange+ii+1],s=area,alpha=0.5,color='red')
                        vv=data1[kpos,rangv-NoteRange+ii+1]
                    for ii in range(NoteRange,2*NoteRange):
                        plt.scatter(data1[kpos,0],data1[kpos,rangv-NoteRange+ii+1],s=area,alpha=0.5,color='blue')
                        cc=data1[kpos,rangv-NoteRange+ii+1]
    print ('---------end--------')
    print (kpath)
    plt.xticks(kvalue,kpath,color='black')
    axes.xaxis.grid()
    if titleTag:
        if setcolor=='none':
            plt.title('Exact band diagram',weight='bold')
        elif setcolor=='all':
            plt.title('Diagrams in one',weight='bold')
        elif n*m<0:
            plt.title('Energy range ('+str(n)+'eV ,'+str(m)+'eV)',weight='bold')
        else:
            plt.title('('+str(n)+' ,'+str(m)+')',weight='bold')
    


plt.rc('font', size=16,weight='bold')# controls default text sizes
plt.rc('axes', titlesize=15)    # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)   # fontsize of the tick labels
plt.rc('ytick', labelsize=15)   # fontsize of the tick labels
plt.rc('legend', fontsize=15)   # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title


fig,axes = plt.subplots(nrows=1, ncols=1,figsize=(10,15))

DataFile='example1.sh'
dataDIR='./BandData/00-'

plt.sca(axes)
taglist = [1]
drawDiagrams(dataDIR,axes,DataFile,taglist,'all',0)
plt.ylabel('Energy (eV)',weight='bold')


fig.tight_layout()
plt.savefig('band_diagrams.png')
plt.show()


