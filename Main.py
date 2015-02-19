'''
Created on Dec 14, 2014

@author: ebadr
'''

import sys
sys.setrecursionlimit(10000)

from ProfileMatrix import PrMt
from MCS import Main
from MCSTree import MainG
from GenerateSeqs import  main

#4. write Frequencies
def Write_Seqs(F,Seqs):
    # fh = None
    fh = open(F, "w", encoding="utf8")
    for s in Seqs:
        for k in s:
            fh.write(str(k))
            fh.write(',')
        fh.write("\n")
    fh.close()

#$$$$$$$$$$$$$$$4444
 
R1=sys.argv[1]
Alpha1=sys.argv[2]
Theta1=sys.argv[3]
GraphW1=sys.argv[4] 
    
def Runall(R1,Alpha1,Theta1,GraphW1):
    print(R1,Alpha1,Theta1,GraphW1)
    R=int(R1)
    Alpha=int(Alpha1)
    Theta=int(Theta1)
    GraphW=str(GraphW1)
    ID=1
    ###################
    print('start at -->')  
    
    print('-----> Profile Matrix part 1')
    M1,G_ESE,ESEs= PrMt(R,GraphW,'ESE1')
    print('-----> Finished Profile Matrix part 1')
    ############################
    print('-----> Profile Matrix part 2')
    M2,G_ESS,ESSs= PrMt(R,GraphW,'ESS1')
    print('-----> Finished Profile Matrix part 2')
    ############################
    print('-----> MCS part 1')
    MH, ID, attributes1 = Main(ID,R,Alpha,M1,G_ESE,GraphW,ESEs,'ESE1')
    print('No. of enhancers = ',str(ID-1))
    ESEID= ID
    print('-----> Finished MCS part 1')
    ############################
    print('-----> MCS part 2')
    mh, ID, attributes2 = Main(ID,R,Alpha,M2,G_ESS,GraphW,ESSs,'ESS1')
    print('No. of Silencers',str(ID-ESEID))
    print('-----> Finished MCS part 2')
    ############################
    #adding WCCs
    MH.extend(mh)
    print('No. of MCSs',len(MH))
    MH1={}
    MCST={}
    for i in range(len(MH)):
        w = MH[i]
        MCST[w.graph['G']]= w.graph['att']
        MH1[w.graph['G']]= w
    
    attributes1.extend(attributes2)
    Write_Seqs(GraphW+'MCSs/attr.txt',attributes1)
    ############################
    print('-----> PrefixTreeGraph')
    MsetGraphs, MsetExons = MainG(MH1,MCST,Theta, GraphW+'MCSs/')
    print(len(MsetGraphs))
    print('-----> FinishedPrefixTreeGraph')
    ##########################
    for j in MsetGraphs:
        L = MsetGraphs[j]
        L2=[int(i) for i in L]
        MsetGraphs[j]= L2 
    ##########################
    print('-----> Generate_MCS_Sequences')
    main(MH1,MsetGraphs, MsetExons,Theta,ESEID,GraphW)
    print('-----> FinishedGenerate_MCS_Sequences')
    return

if __name__== "__main__":
    Runall(R1,Alpha1,Theta1,GraphW1)
#Runall(1000, 50, 20, "/home/ebadr/Writing/py3eg/Tissues/Testis/")
