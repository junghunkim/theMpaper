#!/usr/bin/zsh

g++ -O3 -larmadillo -o myEstimator myEstimator.cpp 

((itr = 0))
for i in DATA/*; do
cp $i myData.txt
./myEstimator 12 3 10 10
cp SIMOUT.mat RESULT/myRESULT_$itr.txt
((itr++))
done