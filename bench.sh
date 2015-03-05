#!/bin/sh
rm -f bench.tmp

for x in {1..10}
do
    time -f "%e" -o output/bench.txt -a $@
done

awk '
BEGIN {max = 0; min = 100000}
    {s+=$0}
    {n+=1}
    {if (max<$0) max=$0}
    {if (min>$0) min=$0}
END {print "avg: " s/n; print "min: " min; print "max: " max}' bench.txt
