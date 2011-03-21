icc -O3 -o myEstimator myEstimator.cpp -I/usr/local/armadillo/ubuntu64/usr/include

R < SIMIN_generator.R --no-save

((itr = 0))
for i in DATA_34v/*; do
cp $i myData.txt
./myEstimator 34 3 50 50 | tee myRESULT_$itr.txt & 
cp SIMOUT.mat RESULT_34v/myRESULT_$itr.txt
((itr++))
done
