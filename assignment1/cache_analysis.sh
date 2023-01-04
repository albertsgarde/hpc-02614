module load studio

rm -rf cache_experiments
mkdir cache_experiments

MS="10 50 100 250 375 500 625 750 875 1000 1125 1250 1500 1750 2000 3000 4000"
TYPES="nat lib mnk mkn nmk nkm kmn knm"

rm -f cache.dat

for t in $TYPES
do
  for m in $MS
  do
    collect -o cache_experiments/${t}_${m}.er -p on -S on -h dcm -h insts -h l2m -h l3m -j off /zhome/ab/8/137185/hpc/hpc-02614/assignment1/matmult_c.gcc $t $m $(($m+10)) $(($m+20))
    er_print -script er_print_script cache_experiments/${t}_${m}.er | grep "matmult_[^r-r][^e-e][^f-f]" >> cache.dat
    truncate -s -1 cache.dat
    echo ",$m" >> cache.dat
  done
done
