# CoSREM

What is CoSREM?
CoSREM (Combinatorial SRE Miner) is a graph mining algorithm 
to discover combinatorial SREs in human exons. 

Requirements:
This package is implemented under Python 3.3. and Centos 6.

It was developped using PyDev.

Python packages to be intalled before use:
networkx
matplotlib
numpy

Using CoSREM:

1- Place the required data(files provided under Data branch) in a specific directory
   Required files are:
	  hexmer-Ei-Order.txt : contains all possible 6-mers with their corresponding scores and ranks.
	  exons+introns+new.csv : contains all unique coding exons for known human gene with thier flanking intronic regions.
	  Exons1.txt : contains the first 50 nucleotides from the utilized exons

2- Create a folder called MCSs in the same directory

3- Run CoSREM from terminal using this command
   python Main.py R Alpha Theta Path_to_directory
   where(R, Alpha, Theta, Path_to_directory) are user-defined parameters.
	 R: the number of 6-mers with the highest ranks to be considere and the number of 6-mers with the lowest ranks as 		    well.
	 Alpha: the number of exons that two 6-mers should at least have in common to be considered one longer k-mer.
	 Theta: the number of exons a set of SREs should at least have in common to be considered combinatorial SRE set.   
	 Path_to_directory: this is a string specifying the path to the input files.

4- After running CoSREM
   The output files:
	  1- All maximal cohesieve subgraphs will be created in the MCSs folder in gml format
	  2- attr.txt is a file that contains the attribues (shared exons IDs) of each maximal cohesive subgraph (MCS)   		     ordered by the MCS IDs.
	  3- FinalSREs.txt includes the final result of all identified combinatorial SREs after the filtering stage.   
 
For further questions please contact 

  Eman Badr  
  
  ebadr@vt.edu




