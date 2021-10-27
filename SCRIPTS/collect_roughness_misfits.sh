elem=Mg #change to current working element

#initiate temporary files:
echo "LAMBDA" > temp_l.txt
echo "ROUGHNESS_SUM" > temp_r.txt
echo "MISFIT" > temp_m.txt

#add data to temp files
for file in $(ls ../DATA/INVERSE_RESULTS/$elem*.txt)
do
    awk 'NR==1{print $1}' $file >> temp_l.txt
    x=`awk 'NR==2{print $1}' $file`
    y=`awk 'NR==3{print $1}' $file`
    echo $x
    echo $y
    rough_sum=`echo $x+$y | bc`
    echo $rough_sum
    echo $rough_sum >> temp_r.txt
    awk 'NR==4{print $1}' $file >> temp_m.txt
done

pr -mts' ' temp_l.txt temp_r.txt > temp_concat.txt
pr -mts' ' temp_concat.txt temp_m.txt > ../DATA/INVERSE_RESULTS/$elem\_tradeoff.txt
#rm temp*