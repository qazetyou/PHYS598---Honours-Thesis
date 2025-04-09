This is the final code submitted as part of my senior thesis.

Package requirement: nauty 2.8.9

Functionalities:
- Find local complements of a given graph.
- Compare the graph to another graph and check if they are in each other's local complementation orbit.
- Find the local complements that satisfy given vertex permutations (check if the vertex permutations are in the graph's automorphism group).
- Count the number of unique graphs that satisfy the given vertex permutations, up to isomorphism.

Description:
The program asks for the number of vertices and edges of the user-given graph. Please note that _the vertex indices in the program starts from 0 instead of 1_ so all the vertex labels obtained from the lattice must be offset by 1. It then asks if the user wants to check if they want to compare this graph to another graph. If the user chooses to do so, they enters the new graph and the program will later confirm whether or not the two graphs are LC-equivalent. The program then asks the number of vertex permutation the user wants to check, and the specific vertex permutations. Finally, the program inquires the user for their prefered depth of the local complementation path. 
After the querries, the program outputs the local complement that satisfies all vertex permutations and the vertex path that leads to it, if any is to be found. If the user has chosen the option to compare graphs, the program also gives the result of this. Finally, the program outputs the number of unique graphs that satisfy all the vertex permutations, up to isomorphism. If the number is greater than one, all unique graphs will be displayed.
