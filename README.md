# kLSVC
A program to perform local search of vertex covers by systematically exploring local search neighborhoods.

# Installation
Simply download the code and call `make` in the directory where the Makefile is located.. You'll find the executable `kLSVC` in the `bin` directory contained in the same folder as the Makefile. The code was confirmed to run on macOS 10.12.2 and Debian GNU/Linux 7.

To simplify analysis there is also a compile option that reduces the output of the program. By using

```make reducedoutput```

the program will output something like

```
k,duration,searchNodes,improvementFound
1,0.036,314,0
3,0.93,1487,1
...
```

which indicates that it took 0.036 seconds to process 314 search nodes at a neighborhood radius of k = 1 and no improvement was found. Afterwards it took 0.93 seconds to process 1487 search nodes at a neighborhood radius of k = 3 and an improvement was found.

#Usage
The program takes a graph and a vertex cover of the graph and tries to improve the given cover using local search. The graph is expected to be a simple text file containing one line for each edge and the vertices being numbered consecutively from 1 to n, where n is the size of the graph.

By not passing a cover the program will compute a greedy vertex cover:

```./kLSVC -graph path/to/graph.txt -out path/to/greedyCover.txt```

In order to perform a local search with a given cover use:

```./kLSVC -graph path/to/graph.txt -cover path/to/cover.txt -kMax 21 -out path/to/improvedCover.txt```

Here -kMax specifies the maximum neighborhood size to check.

Pruning is turned on by passing `-pr2`or `-pr3` (or both):

```./kLSVC -graph path/to/graph.txt -cover path/to/cover.txt -kMax 21 -pr2 -pr3 -out path/to/improvedCover.txt```
