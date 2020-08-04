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
FM-index[[1]](#1) is a commonly used alignment algorithm in state-of-the-art software such as Bowtie. It is a space-efficient data structure that combines the properties of suffix array with the Burrows-Wheeler transform (BWT)~\cite{bwt}. This data structure provides an efficient mechanism to perform substring matching of a pattern $P$ in a long reference sequence $R$. 

The generation of FM-index begins with the computation of BWT of the reference genome `R`, i.e. `BWT(R)`.
First, $R$ is terminated with a unique character: `$`, which is lexicographically the smallest value. Then, all the rotations of the text are obtained and sorted correspondingly. The suffix array can be obtained by considering the characters before `$` in each entry of the rotation list. `BWT(R)` can be formed by extracting and concatenating the last characters of all the entries on the sorted list. \tabref{tab:bwt} demonstrates the derivation of BWT with an example reference genome `R=GCTAT`. The string preceding the `$` sign in the sorted rotations forms the suffix array (`SA`), which indicates the position of each possible suffix in the original string.

## References
<a id="1">[1]</a> 
 P. Ferragina and G. Manzin,
"An Experimental Study of an Opportunistic Index,"
12th Annual ACM-SIAM Symposium on Discrete Algorithms, 2001, pp. 269-278.
