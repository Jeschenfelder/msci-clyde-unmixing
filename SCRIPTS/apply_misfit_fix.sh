for file in $(ls ../DATA/INVERSE_RESULTS/Mg*.asc)
do
    python3 fixing_inverse_misfits.py $file
done