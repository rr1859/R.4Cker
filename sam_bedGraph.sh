#!/bin/bash

directory=$1
cd $directory
for file in *aligned.sam ; do awk < $file '{print $4}' | grep "^ch" | sort | uniq -c | sed -e 's/_\|:\| \|-/\t/g' |  awk '{print $2"\t"$3"\t"$4"\t"$1 }' - >  $(echo $file | sed 's/.sam/.bedGraph/g') ; done
