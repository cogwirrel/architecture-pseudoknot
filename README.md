architecture-pseudoknot
=======================

Advanced Architecture Coursework 1

To get min value of column 4 in awk:

`cat results/ruulsq.out | awk -F " " 'NR == 1 {min = $4} $4 < min {min=$4;minline=$0} END{print minline}'`