#!/bin/bash

dirlist=`find . -maxdepth 1 -type d | grep kl | sed 's/\.\///g'`
echo $dirlist

for dir in $dirlist
do
  echo "Entering $dir"
  cd $dir
  filelist=`ls *.dat`
  for fname in $filelist
  do
    echo "gzip $fname"
    gzip $fname
  done
  echo "Leaving $dir"
  cd ..
done
