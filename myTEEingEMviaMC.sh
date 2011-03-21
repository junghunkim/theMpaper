icc -O3 -o myTEEingEstimator myEstimator.cpp -I/usr/local/armadillo/ubuntu64/usr/include

R < SIMIN_generator.R --no-save

((itr = 0))
for i in DATA_34v/*; do
cp $i myData.txt
./myEstimator 34 3 50 50 | tee myTEEedRESULT_$itr.txt
((itr++))
done
cp myTEEedRESULT_*.txt RESULT_34v/
