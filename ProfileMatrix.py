'''
Created on Dec 14, 2014

@author: ebadr
'''

#$$$$
#!/usr/bin/env python3
#import
import sys
#import math
import networkx as nx
import matplotlib
import numpy as np

matplotlib.use('TkAgg')
try:
    import matplotlib.pyplot as plt
except:
    raise
sys.setrecursionlimit(10000)
#from operator import itemgetter

#Functions
#1. Read the 4096 hexamers and their orders
#Read Excel Sheet where we have kmers and their Enrichement index (I'll make it csv format)
def Read_Data(FileName):
    fh = None
    try:
        fh = open(FileName, encoding="utf8")
        Kmers = []
        Index = {}
        Order = {}
        for lino, line in enumerate(fh, start=1):
            line = line.rstrip()
            if not line:
                continue
            elif "," in line:
                kmer, index, order = line.split(",", 2)
                Kmers.append(kmer)
                order=order[:-1] #activate it if enhancers file
                Index[kmer]=index
                Order[kmer]=order
            else:
                raise KeyError("parsing error on line {0}".format(lino))
        return Kmers, Order
    finally:
            if fh is not None:
                fh.close()

#$$$$$$$$$$$$$$$$$$$$$
def Read_Input_Data(FileName,length):
    fh = None
    try:
        fh = open(FileName, encoding="utf8")
        ID=[]
        Exons1=[]
        Exons2=[]
        Introns1=[]
        Introns2=[]
        tmp=""
        A=""
        p5=200
        p3=200
        rem='N'
        for lino, line in enumerate(fh, start=1):
            line = line.rstrip()
            if not line:
                continue
            else:
                A,tmp = line.split(",", 1)
                intron1=tmp[0:p5]
                intron2=tmp[len(tmp)-p3:len(tmp)]
                tmp=tmp[p5:len(tmp)-p3]
                if len(tmp)>=(length*2):
                    a,b,c,d,e,f= A.split(' ',5)
                    f = 0
                    try:
                        b,c,d,e,f = a.split('_',4)
                        #print(b,c,d,e,f)
                    except ValueError:
                        f = 0
                    
                    f=int(f)+1
                    ID.append(int(f)) 
                    Introns1.append(intron1)
                    Introns2.append(intron2)
                    Exons1.append(tmp[0:length])
                    Exons2.append(tmp[len(tmp)-length:len(tmp)])
                    tmp=""
        return Exons1, Exons2, Introns1, Introns2, ID
    finally:
            if fh is not None:
                fh.close()
#####################
def Read_Input_Data_Tissue(FileName,length):
    fh = None
    try:
        fh = open(FileName, encoding="utf8")
        Exons1=[]
        Exons2=[]
        tmp=""
        A=""
        for lino, line in enumerate(fh, start=1):
            line = line.rstrip()
            if not line:
                continue
            else:
                A,b,c,d,e,f,g,h,tmp = line.split(",", 9)
                if len(tmp)>=(length*2):
                #if len(tmp)>=(length): #tissue
                    Exons1.append(tmp[0:length])
                    Exons2.append(tmp[len(tmp)-length:len(tmp)])
                    tmp=""
        return Exons1, Exons2
    finally:
            if fh is not None:
                fh.close()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#3. Build the whole debruijn Graph and save it in a form so we can upload it again
#Build De bruijn graph. Each node is associated with order and edge is 
def Build_Graph(G,Order,hexa,n):
    # s=nx.number_of_nodes(G)
    if(n==4097):    
        return G
    else:
        bases=('A','C','G','T')
        overlap=hexa[1:len(hexa)]
        for i in range(4):
            new_hexa= overlap+bases[i]
            if new_hexa in G:
                G. add_edge(hexa,new_hexa)
            else:
                G.add_node(new_hexa, order = Order[new_hexa])
                #G.add_node(new_hexa)
                G.add_edge(hexa,new_hexa)
                n=n+1;
                Build_Graph(G,Order,new_hexa,n)
    return G
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#4. Build P matrix
def Build_Matrix(SREs,Exons):
    S=len(SREs)
    #P = np.zeros((S,2*len(Exons1)), dtype=np.int)
    P = np.zeros((S,len(Exons)), dtype=np.int)
    #print(P.shape)
    # Enhancers and silencers
    for i in range(S):
        K=SREs[i]
        for j in range(0,len(Exons)):
            E=Exons[j]
            if (K in E):
                P[i,j]=1
    return P

#$$$$$$$$$$$$$$$4444
#4. write Frequencies
def Write_Seqs(F,Seqs):
    # fh = None
    fh = open(F, "w", encoding="utf8")
    for s in Seqs:
        fh.write(str(s))
        fh.write("\n")
    fh.close()

#$$$$$$$$$$$$$$$4444



def PrMt(R,GraphW,flag):    
#Variables
    length=50 #nucleotides
    SREs=[]
    filelist= [f for f in os.listdir(GraphW+'MCSs') if f.endswith(".gml")]
    for f in filelist:
        os.remove(GraphW+'MCSs/'+f)
    FileName= GraphW +"hexmer-Ei-Order.csv"
    File= GraphW +"exons+introns+new.csv" #tissue
    #File = GraphW +"Testis.csv"
    Ex = GraphW+"Exons1.txt"
#Calling Functions

    #print("Calling Functions")
    Kmers, Order = Read_Data(FileName)
    Exons1,Exons2, Introns1, Introns2,ID= Read_Input_Data(File,length) #tissue
    #Exons1,Exons2= Read_Input_Data_Tissue(File,length)
    #Write_Seqs(Ex,Exons1) #tissue
    if flag== "ESE1":
        for i in range(R):
            SREs.append(Kmers[i])
    elif flag== "ESS1":
        for i in reversed(range(R)):
            SREs.append(Kmers[len(Kmers)-i-1])
    G=nx.DiGraph()
    G.add_node('AAAAAA', order = Order['AAAAAA'])
    G=Build_Graph(G,Order,'AAAAAA',0) 
    P=Build_Matrix(SREs,Exons1)
    Gu=nx.subgraph(G,SREs)
    return P,Gu,SREs

#PrMt(400,'/home/ebadr/Writing/py3eg/','ESE1')
