lru=(4 8 16 24 32 40 48 56 64 128)
lfu=(128)

echo "BEGIN"
for i in ${lfu[@]}; do
    for j in ${lru[@]}; do
        echo "LRU $j, LFU $i"
        ./out/GenomeClustering res/p3_clean_C-148-2-Caecum_S128_sorted.fa clusters/clusters_mufand_${j}_${i}.uc 0.03 -r $j -f $i
    done
done
echo "DONE"
