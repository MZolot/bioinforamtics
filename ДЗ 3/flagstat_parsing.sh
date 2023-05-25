#!/usr/bin/bash

l=$(head -8 results/flagstat_result.txt | tail -1)

symbol=d
p1=$(expr index "$l" $symbol)
p2=$[p1+2]
numb=${l:p2:1}

if [ "$numb" -eq 7 ]
then 
echo "OK"
else echo "NOT OK"
fi