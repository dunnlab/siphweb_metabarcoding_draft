#!/bin/bash

splitfiles=($(ls "splitout"))
for i in ${splitfiles[@]}; do

echo ~/anaconda2/bin/sap --project splitSAPs/$i --email alejandro.damianserrano@yale.edu --minsignificance 0.3 -s 0.1 --minidentity 0.5 -n 10 --ppcutoff 85 --svg splitout/$i

done
