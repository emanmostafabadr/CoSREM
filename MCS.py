'''
Created on Dec 14, 2014

@author: ebadr
'''
'''
Created on Sep 4, 2014

@author: ebadr
'''

import sys
import os
import copy
#import math
import networkx as nx
#import matplotlib
import csv
#import numpy as np
#import Bio as bi
#from Bio import trie

#matplotlib.use('TkAgg')
#try:
#    import matplotlib.pyplot as plt
#except:
#    raise
sys.setrecursionlimit(10000)
#from operator import itemgetter

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

def Read_Input_Data(FileName):
    fh = None
    try:
        fh = open(FileName, encoding="utf8")
        SREs=[]
        for lino, line in enumerate(fh, start=1):
            line = line.rstrip()
            if not line:
                continue
            else:
                SREs.append(line)
        return SREs
    finally:
            if fh is not None:
                fh.close()
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def Sum(X):
    count=0
    for i in range(len(X)):
        count += X[i]
    return count
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def Cohesive_Single_Genes(SREs,M,Alpha):
    V=[]
    for i in range(len(SREs)):
        #s= Sum(M[i,:])
        s= Sum(M[i])
        #print(s)
        #s= np.sum(M[i,:])
        #print(s)
        if(s>=Alpha):
            #   print(s)
            V.append(SREs[i])
    V.sort()
    return V
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#common attributes between the nodes of a subgraph
def Attributes(sre,row):
    att=[]
    for i in range(len(row)):
        if(row[i]==1):
            att.append(i)
    return att
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def check(G,MH):
    for elem in MH:
        e= nx.nodes(elem)
        n=nx.nodes(G)
        if set(n) == set(e):
            # attributes
            if ( set(elem.graph['att']) == set(G.graph['att'])):
                return True
    return False

#$$$$$$$$$$$$$$$$$$$$$$$$$$$44
def Get_Neighbours(G,Gu):
    Ne=[]
    tmp=[]
    n = nx.nodes(G)
    for e in n:
        s= nx.all_neighbors(Gu, e)
        #s= Gu.neighbors(e)
        tmp= set(s) | set(tmp)
    Ne = list(tmp-set(n)) 
    return Ne
#$$$$$$$$$$$$$$$$$$$$$$$$$$$44
def genMCS(G,Alpha,MH,S,M,Gu,Order,SREtype):
    if G in S:
        return
    else:
        S.append(G)
        #print("S",S)
    f = check(G,MH)
    if(f):
        return
    maxFlag=True
    #Get the neighbors of the graph
    Ne =  Get_Neighbours(G,Gu)
    Ne.sort()
    for elem in Ne:
        nodes= nx.nodes(G)
        nodes.append(elem)
        NewG=nx.DiGraph(Gu.subgraph(nodes))
        #Get the attributes of the new graph by intersection of the old graph and only the one node that was added
        #NewG.graph['att']=list(set(G.graph['att'])& set(Attributes(int(Order[elem])-1-3696,M[int(Order[elem])-1-3696,:])))
        if SREtype == 'ESE1' or SREtype == 'ESE2' or SREtype == 'ESER':
            #NewG.graph['att']=list(set(G.graph['att'])& set(Attributes(int(Order[elem])-1,M[int(Order[elem])-1,:])))
            NewG.graph['att']=list(set(G.graph['att'])& set(Attributes(int(Order[elem])-1,M[int(Order[elem])-1])))
        else:
            #NewG.graph['att']=list(set(G.graph['att'])& set(Attributes(int(Order[elem])-1-3696,M[int(Order[elem])-1-3696,:])))
            NewG.graph['att']=list(set(G.graph['att'])& set(Attributes(int(Order[elem])-1-3696,M[int(Order[elem])-1-3696])))
        # check if the new graph is cohesive
    #    print("NewG",len(NewG.graph['att']))
        if (len(NewG.graph['att'])>=Alpha):
            maxFlag=False
            genMCS(NewG,Alpha,MH,S,M,Gu,Order,SREtype)
    if maxFlag == True:
        #print("G",G.nodes())
        # print("-->",len(G))
        MH.append(G)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def Generate_MCS(R,V,M,Alpha,Gu,Order,SREtype):
    S=[]
    MH=[]  
    #for sre in V:
    for i in range(len(V)):
        sre=V[i]
        r=[]
        r.append(sre)
        G = Gu.subgraph(r).copy()
        if SREtype == 'ESE1' or SREtype == 'ESE2' or SREtype == 'ESER':
            #G.graph['att']=Attributes(int(Order[sre])-1,M[int(Order[sre])-1,:]) #enhancers
            G.graph['att']=Attributes(int(Order[sre])-1,M[int(Order[sre])-1]) #enhancers
        else: 
            #G.graph['att']=Attributes(int(Order[sre])-1-3696,M[int(Order[sre])-1-3696,:]) #silencers 
            index=4096-R
            #G.graph['att']=Attributes(int(Order[sre])-1-3696,M[int(Order[sre])-1-3696]) #silencers
            G.graph['att']=Attributes(int(Order[sre])-1-index,M[int(Order[sre])-1-index]) #silencers 
        genMCS(G,Alpha,MH,S,M,Gu,Order,SREtype)
    return MH
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4   


#4. write Frequencies
def Write_Seqs(F,Seqs):
    # fh = None
    fh = open(F+'attr.txt', "w", encoding="utf8")
    for s in Seqs:
        for k in s:
            fh.write(str(k))
            fh.write(',')
        fh.write("\n")
    fh.close()

#$$$$$$$$$$$$$$$4444



#Alpha = 1000


def Main(ss,R,Alpha, M, Gu, Path, SREs, SREtype):
    
    #Variables
    FileName= Path+"hexmer-Ei-Order.csv"    
    GraphW= Path+"MCSs"  
    ID=ss
###################################################################33

    #Calling Functions
    Kmers, Order = Read_Data(FileName)
    # Start pruning the nodes that doesnt have the minimum attributes
    V = Cohesive_Single_Genes(SREs,M,Alpha)
    #print('V',len(V))
    MH = Generate_MCS(R,V,M,Alpha,Gu,Order,SREtype)
    #Write the MCSs
    GraphW= GraphW+'/MCS'
    attributes=[]
    for i in range(len(MH)):
        w = MH[i]
        w.graph['G']=ss
        #print(w.nodes(),w.graph['G'],len(w.graph['att']))
        attributes.append(list(w.graph['att']))
        ss=ss+1
    
    MH1 = copy.deepcopy(MH)
    
    for i in range(len(MH)):
        w = MH[i]
        del w.graph['att']
        nx.write_gml(w,GraphW+str(ID)+'.gml')
        ID=ID+1
    
    return  MH1, ss ,attributes 
