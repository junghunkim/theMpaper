#!/usr/bin/zsh
for i in param_init_dir/*; do
cp $i ./Param_Init.R
R < datagen.r --no-save
cp myoutfile.txt ${i:s/.R/_output.txt}
cp myoutfile.txt myData.txt
./param_estimator 4 3 10 10
cp SIMOUT.mat ${i:s/.R/_filtered.mat}
cp SIMOUT.mat SIMIN.mat 
done
