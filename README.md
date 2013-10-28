## 3-node Network Motifs ##
> Shane Neph


Overview
=========
A set of 'network motif' programs.  

Briefly, a network motif is a particular configuration of directed edges between any N-node subnetwork/circuit.  For the three-node subnetwork case, there are 13 distinct configurations possible between connected nodes.  One can break down any directed graph into a collection of these distinct configurations.  A natural question to ask is whether you find more (or fewer) instances of a network motif in a graph than you would expect by chance, given a directed graph of the same number of nodes and edges.

Research has shown the theoretical and measured importance of each of these simpler three-node circuits, and that directed graphs found from a variety of disciplines often share signatures over the relative enrichment/depletion of these simple circuits.  

See, for example, Network Motifs: Simple Building Blocks of Complex Networks, Science, 298:824-827 (2002)  


Programs
=========
find_3node_motifs  
Quickly enumerate all 3-node network motifs in a directed graph.  

motif3_network_changes  
What do all three-node network motifs look like in a different graph?  This is useful for comparing two graphs.  Say you have biological regulatory networks in graph form over two cell type samples.  How do all 3-node feedforward loops found in one cell type's graph distribute over all 13 three-node network motifs in the second cell type's graph?

Build
======
// requires g++ version 4.7 or newer  
make -C src/

How-To
=======
find_3node_motifs \<input-graph\> \> output.results  
  \<input-graph\> is a file with rows of the form:  
A   B  
  where a tab separates the node labels and A->B in your graph  


motif3_network_changes \<target-network-file\> \<reference-network-file\> \> output.mtx  
  both graphs should contain rows of the form:  
A   B  
  where a tab separates the node labels and A->B in your graph.

  Determine how every 3-node circuit in \<reference-network-file\> is configured over the same nodes in \<target-network-file\>.  'not-found' is an additional category when the 3 nodes are not connected in \<target-network-file\>.
