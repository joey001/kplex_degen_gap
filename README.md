# Fast Maximum $k$-Plex Algorithms Parameterized by Small Degeneracy Gaps

The work extends the conference paper "A Fast Maximum $k$-Plex Algorithm Parameterized by the Degeneracy Gap", presented in the 32nd International Joint Conference on Artificial Intelligence (IJCAI'23) in Macao on 19th August 2023.

## Files

* code - codes of our algorithms and kPlexS and KpLeX.
* data - example graph instances
* supplement - supplementary materials

## Compile the code

```sh
$ make
```
It generates the executable (kPlexS, KpLeX, Maple and its variants).

## Run the code

```sh
$ ./EXE {path_to_binary_compressed_graph} {k_value}
```

An example of computing the exact maximum 15-plex for the dataset soc-slashdot is as follows
```sh
$ ./Maple data/soc-slashdot.bin 15
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

Note that:

The lower bound of both kPlexS and KpLex are simply set to at least 2k-2 to screen out trivial cases.

And the input of kPlexS and KpLeX are modified to accept binary encoded instances.