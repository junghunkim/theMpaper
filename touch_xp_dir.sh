#!/usr/bin/zsh
((num = 1))
mkdir ./xp_dir_$num
mkdir ./xp_dir_$num/param_init_dir

cp ./run_mpaper.sh ./xp_dir_$num/run_mpaper.sh
cp ./datagen.r ./xp_dir_$num/
cp ./simulator.r ./xp_dir_$num/

for ((i=0;i<10;i++)); do
cp ./Param_Init.R ./xp_dir_$num/param_init_dir/Param_Init_$i.R
done
