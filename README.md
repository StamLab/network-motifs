## 3-node Network Motifs ##
> Shane Neph


Overview
=========
A set of 'network motif' programs.  

Briefly, a network motif is a particular configuration of directed edges between any N-node subnetwork/circuit.  For the three-node subnetwork case, there are 13 distinct configurations possible between connected nodes.  One can break down any directed graph into a collection of these distinct configurations.  A natural question to ask is whether you find more (or fewer) instances of a network motif in a graph than you would expect by chance, given a directed graph of the same number of nodes and edges drawn at random (a Erdős–Rényi graph).  A more appropriate null model might further require that the nodes have the same number of incoming/outgoing edges as in your directed graph.

Research has shown the theoretical and measured importance of each of these simpler three-node circuits, and that directed graphs found from a variety of disciplines often share signatures over the relative enrichment/depletion of these simple circuits.  

See, for example, Network Motifs: Simple Building Blocks of Complex Networks, Science, 298:824-827 (2002)  


Programs
=========
_find_3node_motifs_  
Quickly enumerate all 3-node network motifs in a directed graph, with connection-type information. 

_motif3_network_changes_  
What do all three-node network motifs look like in a different graph?  This is useful for comparing two graphs.  Say you have biological regulatory networks in graph form over two cell type samples.  How do all 3-node feedforward loops found in one cell type's graph distribute over all 13 three-node network motifs in the second cell type's graph?

Build
======
// requires g++ version 4.7 or newer  
make -C src/

How-To
=======
_find_3node_motifs_ [input-graph] \> output.results  
  [input-graph] is a file with rows of the form:  
A   B  
  where a tab separates the node labels and A->B in your graph  


_motif3_network_changes_ [target-network-file] [reference-network-file] depends upon outputs from _find_3node_motifs_  

find_3node_motifs graph-A \> output.graphA  
find_3node_motifs graph-B \> output.graphB  

motif3_network_changes output.graphA output.graphB \> output.mtx  

  Determines how every 3-node circuit in output.graphB is configured over the same nodes in output.graphA.  'No-Match' is an additional category when the 3 nodes are not connected in [target-network-file].  There is also a 'Matched-Variant' column which shows the number of times a circuit is the same between networks, but the arrows between the 3 nodes have changed directions.
