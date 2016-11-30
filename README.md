# kLSVC
A program to perform local search of vertex covers by systematically exploring local search neighborhoods.

# Installation
Simply download the code and call make. The code was confirmed to run on macOS and Debian GNU/Linux 7.

#Usage
The program takes a graph and a vertex cover of the graph and tries to improve the given cover using local search. The graph is expected to be a simple text file containing one line for each edge and the vertices being numbered consecutively from 1 to n, where n is the size of the graph.

By not passing a cover the program will compute a greedy vertex cover:

```./kLSVC -graph path/to/graph.txt -out path/to/greedyCover.txt```

In order to perform a local search with a given cover use:

```./kLSVC -graph path/to/graph.txt -cover path/to/cover.txt -kMax 21 -out path/to/improvedCover.txt```

Here -kMax specifies the maximum neighborhood size to check. 

Pruning is turned on by passing `-pr2`or `-pr3` (or both):

```./kLSVC -graph path/to/graph.txt -cover path/to/cover.txt -kMax 21 -pr2 -pr3 -out path/to/improvedCover.txt```
