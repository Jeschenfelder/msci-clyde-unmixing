for file in $(ls ../DATA/INVERSE_RESULTS/)
do
    python3 fixing_inverse_files.py $file
done
