## 3-node Network Motifs ##
> Shane Neph


Overview
=========
A set of 'network motif' programs.  

Network Motifs: Simple Building Blocks of Complex Networks, Science, 298:824-827 (2002)  

find_3node_motifs  
Quickly enumerate all 3-node network motifs in a directed graph.  

motif3_network_changes  
Determine how two graphs differ in all their 3 node network motifs  

Build
======
// requires g++ version 4.7 or newer  
make -C src/

How-To
=======
find_3node_motifs \<input-graph\>  
  \<input-graph\> is a file with rows of the form:  
A   B  
  where a tab separates the node labels and A->B in your graph  


motif3_network_changes \<target-network-file\> \<reference-network-file\>  
  both graphs should contain rows of the form:
A   B
  where a tab separates the node labels and A->B in your graph.

  How do those 3-node network motifs in \<reference-network-file\> map onto the same nodes in \<target-network-file\>?
