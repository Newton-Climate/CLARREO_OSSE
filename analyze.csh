#!/bin/tcsh -f
 
set file = $1

rm -rf before2after.txt
rm -rf after.txt
rm -rf list_o_pnts.txt

grep "BEFORE" $file | sort | sed "s/BEFORE/AFTER/" > before2after.txt
grep "AFTER" $file | sort > after.txt

diff -wbi before2after.txt after.txt | egrep -v "," | egrep -v "d" > list_o_pnts.txt
