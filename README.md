# FPGA-alike software for short read alignment using FM-index 

In the repository `EMBothWay`, DNA alignment (exact-matching) using FM-index is accelerated on Xilinx VU9P using MaxJ. This software aims at generating the corresponding FM-index, and providing FPGA-alike software for testing and debugging purposes.

## Index Generation
```
# cd SW/fm-sa-build
# .fm-sa-build <reference_fasta> <index_name_you_like> true
```
8 index files are generated. <index_name_you_like>.1.ral to <index_name_you_like>.8.ral

## Exact-Match Alignment
```
# cd ../SW/fm-align
# .fm-sa-build ../fm-sa-build/<index_name_you_like> <read>.fq
```
Using the 8 index files, short read alignment can be performed.

## Principle
### FM-Index
FM-index[[1]](#1) is a commonly used alignment algorithm in state-of-the-art software such as Bowtie. It is a space-efficient data structure that combines the properties of suffix array with the Burrows-Wheeler transform (BWT)[[2]](#2). This data structure provides an efficient mechanism to perform substring matching of a pattern `P` in a long reference sequence `R`. 

The generation of FM-index begins with the computation of BWT of the reference genome `R`, i.e. `BWT(R)`.
First, `R` is terminated with a unique character: `$`, which is lexicographically the smallest value. Then, all the rotations of the text are obtained and sorted correspondingly. The suffix array can be obtained by considering the characters before `$` in each entry of the rotation list. `BWT(R)` can be formed by extracting and concatenating the last characters of all the entries on the sorted list. The following table demonstrates the derivation of BWT with an example reference genome `R=GCTAT`. The string preceding the `$` sign in the sorted rotations forms the suffix array (`SA`), which indicates the position of each possible suffix in the original string.

index *n* | SA | Sorted Rotations
------| ---| -----------
0 | 5 |  $GCTAT
1 | 3 |  **AT**$GCT
2 | 1 |  **CTAT**$G
3 | 0 |  **GCTAT**$
4 | 4 |  **T**$GCTA
5 | 2 |  **TAT**$GC

After having the suffix array, `BWT(R)` is sorted to form the `i` and `c` tables. For each element `x` of the alphabet of `R`, `i(x)` is defined as the index of its first occurrence in `sorted-BWT(R)`. For each index `n` in `BWT(R)` and for each character `x` in the alphabet, `c(n, x)` stores the number of occurrences of `x` in `BWT(R)` in the range `[0, n-1]`. The following table illustrates the `i(x)` and `c(n, x)` tables for the sequence `R`.

`c(n, x)`
index *n* | A | C | G | T
----------| - | - | - | -
0 | 0 | 0 | 0 | 0
1 | 0 | 0 | 0 | 1
2 | 0 | 0 | 0 | 2
3 | 0 | 0 | 1 | 2
4 | 0 | 0 | 1 | 2
5 | 1 | 0 | 1 | 2
6 | 1 | 1 | 1 | 2

`i(x)`
A | C | G | T
-- | -- | -- | --
1 | 2 | 4 | 5

## References
<a id="1">[1]</a> 
 P. Ferragina and G. Manzin,
"An Experimental Study of an Opportunistic Index,"
12th Annual ACM-SIAM Symposium on Discrete Algorithms, 2001, pp. 269-278.

<a id="2">[2]</a> 
M. Burrows and D. Wheeler,
"Block-sorting Lossless Data Compression Algo-rithm,"
Digital Equipment Corporation, Tech. Rep., 1994.
