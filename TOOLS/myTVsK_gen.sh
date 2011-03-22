#!/usr/bin/zsh

((maxFILEnumbers=10000))
mkdir TEMP
for ((i=0;i<$maxFILEnumbers;i++)); do
echo set.seed\($i\) > mySEEDer.R
echo filename=\"myTVsK_data_$i\" > myFILEer.R
R CMD BATCH ./datagen.r
mv myTVsK_data_$i TEMP/
echo On $i out of $maxFILEnumbers
done
