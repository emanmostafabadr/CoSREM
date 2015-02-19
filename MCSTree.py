'''
Created on Dec 14, 2014

@author: ebadr
'''
'''
Created on Oct 10, 2014

@author: ebadr
'''

'''
Created on Oct 1, 2014

@author: ebadr
'''
import sys
#import math
import networkx as nx
#import matplotlib
#import numpy as np
import operator
#from operator import itemgetter

#try:
#    import matplotlib.pyplot as plt
#except:
#    raise
sys.setrecursionlimit(10000)
##########################################################
def Readseqs(GraphW):
    fh = None
    #ESE1
    try:
        fh = open(GraphW+'attr.txt', encoding="utf8")
        seqs=[]
        for lino, line in enumerate(fh, start=1):
            line = line.rstrip()
            if not line:
                continue
            else:
                L=[]
                L = line.split(",")
                A = [x for x in L if x != '']
                L2=[int(i) for i in A]
                seqs.append(L2)
        #return seqs
    finally:
            if fh is not None:
                fh.close()
    #ESE2
    try:
        fh = open(GraphW1+'attr.txt', encoding="utf8")
        #seqs=[]
        for lino, line in enumerate(fh, start=1):
            line = line.rstrip()
            if not line:
                continue
            else:
                L=[]
                L = line.split(",")
                A = [x for x in L if x != '']
                L2=[int(i) for i in A]
                seqs.append(L2)
        return seqs
    finally:
            if fh is not None:
                fh.close()            
#######################################################
def Find_Graphs2(atr,prev_mcs,MCSt):
    mcs=[]
    exons=MCSt[atr]
    for e in exons:
        if (e in prev_mcs):
            mcs.append(e)
    return mcs
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def check(name,mcs,G):
    curr = name.split('_') # current path
    for node in G:
        prev=node.split('_')
        if set(curr)<set(prev): # these nodes are subset of a previously discovered nodes
            if set(mcs) == set(G.node[node]['Exons']):
                return True       
    return False      
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4   
def add_one_tree(MH,G,i,attri,MCSt,name,Theta):
    prev_name=name
    for k in range(i+1,len(attri)):
        prev_mcs = list(G.node[prev_name]['Exons'])
        mcs = Find_Graphs2(attri[k],prev_mcs,MCSt)
        if (len(mcs)>= Theta):
            name=prev_name+'_'+str(attri[k])
            f = check(name,mcs,G)
            if(f):
                return
            #print(name)
            #w= MH[attri[k]]
            #print(w.nodes())
            G.add_node(name,graph = attri[k],Exons=mcs)
            G.add_edge(prev_name,name)
            add_one_tree(MH,G,k,attri,MCSt,name,Theta)
    return
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4   
def Find_Graphs(atr,MCSt):
    #mcs=[]
    #for Graph in MCSt:
        #if atr in MCSt[Graph]:
            #mcs.append(Graph)
    mcs=MCSt[atr]
    return mcs
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4   
def Build_Prefix_Tree(MH,attri,MCSt,Theta):
    #build the directed graph
    G=nx.DiGraph()
    G.add_node('#',graph=-1, Exons=[])
    for i in range(len(attri)):
        mcs= Find_Graphs(attri[i],MCSt) # find attributes associated with graphs
        if (len(mcs)> Theta): # number of shared exons
            #print('first->'+str(attri[i]))
            name=str(attri[i])
            G.add_node(name,graph = attri[i],Exons=mcs)
            G.add_edge('#',name)
            add_one_tree(MH,G,i,attri,MCSt,name,Theta)
    return G
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4   
def Generate_Prefix(MH,MCSt,Theta):
    # Get all attributes
    #attri=[]
    #for Graph in MCSt:   
        #attri=list(set(MCSt[Graph])| set(attri))
        
    attri=MCSt.keys()    
    # sort the unique attributes
    attri = sorted(attri)
    #start building the prefix tree
    G = Build_Prefix_Tree(MH,attri,MCSt,Theta)
   # nx.write_dot(G,GraphW+'test.dot')
   
    #nx.draw_spring(G,node_size=1000,node_color='r',font_size=10)
    #plt.show()  
    #nx.write_gml(G,GraphW+'graph.gml')
    return G 
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4 
def checked(MsetGraphs,MsetExons,L,M):
    for i in MsetGraphs:
        if set(L) < set(MsetGraphs[i]):
            if set(M)== set(MsetExons[i]):
                return False
    return True
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4 
def Generate_MCSs (G):
    MsetGraphs={}
    MsetExons={}
    k=1
    length=nx.single_source_shortest_path_length(G,'#')
    sorted_x = sorted(length.items(), key=operator.itemgetter(1))
    for i in range(len(sorted_x)-1,0,-1):#elem in length:
        elem,tall=sorted_x[i]
        if tall >= beta:
            L = elem.split('_')
            M = list(G.node[elem]['Exons'])
          #  f = checked(MsetGraphs,MsetExons,L,M)
            f=True
            if (f):
                MsetGraphs[k] = L
                MsetExons[k] = M
                k+=1
    #for i in MsetGraphs:
        #print('number of MCSs:', length[elem],'|', 'number of Exons:',len(M) ,'|', 'MCSs are:', L,'|', 'Exons are:',M)
        #print(i,'|',len(MsetGraphs[i]),'|',len(MsetExons[i]) ,'|', list(MsetGraphs[i]))#,'|',list(MsetExons[i]))
    
    return MsetGraphs,MsetExons

def Write_output(F,MsetGraphs, MsetExons):
    fh = open(F+'MCSsets.txt', "w", encoding="utf8")
    # write MCSset ID | subgraph IDs | Exons IDs 
    for i in range(0,len(MsetGraphs)):
        fh.write(str(i+1))
        fh.write('|')
        L = list(MsetGraphs[i+1])
        M = list(MsetExons[i+1])
        for k in L:
            fh.write(str(k))
            fh.write(',')
        fh.write('|')
        for k in M:
            fh.write(str(k))
            fh.write(',')
        fh.write("\n")
    return

#######################################################
#GraphW="/home/ebadr/Writing/py3eg/WCCs/1000-100/"
GraphW1="/home/ebadr/Writing/py3eg/WCCs/1000-100/"

global beta 
beta = 2  # of graphs in one MCS set

#global theta 
#theta = 100  # of shared exons

def MainG(MH,MCSt,Theta,GraphW):
    #Calling Functions
    print("Calling Functions")   
    
    #Generate MCS tree
    G = Generate_Prefix(MH,MCSt,Theta)
    # Get the depths in order to get only the paths with length more than theta
    MsetGraphs, MsetExons = Generate_MCSs(G)
    Write_output(GraphW, MsetGraphs, MsetExons)
    return MsetGraphs, MsetExons
