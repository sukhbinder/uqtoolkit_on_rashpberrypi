#!/bin/bash

tar xzf qdata.tgz

kind="1 2 3 4 5 6"

touch err.log; rm -rf err.log; touch err.log
for i in $kind
do

  paste newq${i}.dat oldq${i}.dat > tmp.dat
  awk '{print $1-$3" "$2-$4}' tmp.dat > kind${i}.dat
  awk 'function abs(value) { return (value<0?-value:value);} \
       {if ( abs($1) > 1.e-12 ) {print "Error in quad points for kind="j}}' j=$i kind${i}.dat >> err.log
  awk 'function abs(value) { return (value<0?-value:value);} \
       {if ( abs($2) > 1.e-12 ) {print "Error in quad weights for kind="j}}' j=$i kind${i}.dat >> err.log

done

# vandermonde
i=10
paste newq${i}.dat oldq${i}.dat > tmp.dat
awk '{print $1-$3" "$2-$4}' tmp.dat > kind${i}.dat
awk 'function abs(value) { return (value<0?-value:value);} \
     {if ( abs($1) > 1.e-12 ) {print "Error in quad points for kind="j}}' j=$i kind${i}.dat >> err.log
awk 'function abs(value) { return (value<0?-value:value);} \
     {if ( abs($2) > 1.e-12 ) {print "Error in quad weights for kind="j}}' j=$i kind${i}.dat >> err.log

# generic recursive
i=11
paste newq${i}.dat oldq${i}.dat > tmp.dat
awk '{print $1-$3" "$2-$4}' tmp.dat > kind${i}.dat
awk 'function abs(value) { return (value<0?-value:value);} \
     {if ( abs($1) > 1.e-12 ) {print "Error in quad points for kind="j}}' j=$i kind${i}.dat >> err.log
awk 'function abs(value) { return (value<0?-value:value);} \
     {if ( abs($2) > 1.e-12 ) {print "Error in quad weights for kind="j}}' j=$i kind${i}.dat >> err.log

nerr=`wc -l err.log | awk '{print $1}'`

if [ "$nerr" -gt 0 ]; then
  echo " - Found errors in gq.cpp -> check kind*.dat files"
else
  echo " - Quadrature tests passed"
  /bin/rm -rf tmp*.dat kind*.dat newq*.dat
  /bin/rm -rf err.log
  tar czf qdata.tgz oldq*.dat
  /bin/rm -rf oldq*.dat
fi
 
