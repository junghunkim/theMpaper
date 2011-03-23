g++ -O3 -o myEstimator myTEEingEstimator.cpp -I/usr/local/armadillo/ubuntu64/usr/include

R < SIMIN_generator.R --no-save

((itr = 0))
for i in DATA_34v/*; do
cp $i myData.txt
./myEstimator 70 3 50 1 | tee myTEEedRESULT_$itr.txt
((itr++))
done
cp myTEEedRESULT_*.txt RESULT_34v/
