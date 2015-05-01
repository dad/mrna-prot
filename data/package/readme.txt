Datasets accompanying:
   "Accounting for experimental noise reveals that transcriptional control,
   amplified by post-transcriptional regulation, sets steady-state protein levels in yeast,"
   Gábor Csárdi, Alexander Franks, David S. Choi, Eduardo M. Airoldi, and D. Allan Drummond, PLoS Genetics, 2015

Each dataset is tab-delimited and individually annotated.

To read datasets into R (http://www.r-project.org):

x <- read.table("scer-mrna-protein-absolute-estimate.txt", header=T, sep='\t', quote='')

> nrow(x)
[1] 5887
> x[1:10,]
     orf      gene     mrna        prot sd.mrna sd.prot n.mrna n.prot
1  Q0045      COX1   0.1835      18.436   0.333   0.371      3      2
2  Q0050       AI1   0.0967       4.667   0.414   0.480      2      1
3  Q0055       AI2   0.0655       3.942   0.462   0.533      2      1
4  Q0060       AI3   0.0669       2.385   0.443   0.541      2      1
5  Q0065       AI4   0.0666       2.183   0.394   0.528      2      0
6  Q0070 AI5_ALPHA   0.0434       1.021   0.472   0.610      2      0
7  Q0075  AI5_BETA   0.0707       2.114   0.404   0.474      2      0
8  Q0080      ATP8   0.1143       5.350   0.393   0.464      2      0
9  Q0085      ATP6   0.1721      17.272   0.325   0.451      3      1
10 Q0105       COB   0.1559      11.131   0.360   0.467      3      1
