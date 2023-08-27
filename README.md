# A Fast Maximum k-Plex Algorithm Parameterized by the Degeneracy Gap

Our code is based on Chang et al.'s work on Proc. VLDB Endow. 16(2), (2022). 
Thanks for their selfless disclosure of their source code.

Note that: 

In the implementation, there is no need to explicitly build complement graphs , as we use the adjacency matrix to store subgraphs.

The lower bound of both kPlexS and KpLex are simply set to at least 2k-2 to screen out trivial cases.

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "KPLEX", which actually corresponds to the our algorithm for optimization problem (Ours).

## Run the code

```sh
$ ./KPLEX {path_to_binary_compressed_graph} {k_value}
```

An example of computing the exact maximum 15-plex for the dataset soc-slashdot is as follows
```sh
$ ./KPLEX data/soc-slashdot.bin 15
```

The solution is in output file "kplexes.txt".

## Data format
We adopt the time-efficient binary format rather than the text format.  Several demonstration graphs are given in "data" folder.

Transforming [network-repo graphs](http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz) from text format to binary format is as follows:
```sh
$ g++ ./toBin.cpp -o toBin
$ mv ./data/soc-slashdot ./data/soc-slashdot.clq
$ ./toBin data/soc-slashdot.clq
$ ls data/soc-slashdot.bin
```