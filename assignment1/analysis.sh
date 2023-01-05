BIN=$(realpath $1)
OUT_PATH=$(realpath $2)
TYPES=$3
MS=$4
BLOCK_SIZES=$5
EXPERIMENT_DIR=$(realpath $6)

rm -rf $OUT_PATH
mkdir $OUT_PATH

module load studio

rm -rf $EXPERIMENT_DIR
mkdir $EXPERIMENT_DIR

COLLECT_CMD="collect -p on -S on -h dch -h dcm -h l2h -h l2m -h l3h -h l3m -j off"


for t in $TYPES
do
  for m in $MS
  do
    $BIN $t $m $(($m+10)) $(($m+20)) >> $OUT_PATH/perf.dat
    truncate -s -1 $OUT_PATH/perf.dat
    echo " $m" >> $OUT_PATH/perf.dat

    $COLLECT_CMD -o $EXPERIMENT_DIR/${t}_${m}.er $BIN $t $m $(($m+10)) $(($m+20))
    er_print -script er_print_script $EXPERIMENT_DIR/${t}_${m}.er | grep "matmult_[^r-r][^e-e][^f-f]" >> $OUT_PATH/cache.dat
    truncate -s -1 $OUT_PATH/cache.dat
    echo ",$t,$m,1" >> $OUT_PATH/cache.dat
  done
done

for bs in $BLOCK_SIZES
do
  for m in $MS
  do
    $BIN blk $m $(($m+10)) $(($m+20)) $bs >> $OUT_PATH/perf.dat
    truncate -s -1 $OUT_PATH/perf.dat
    echo " $m" >> $OUT_PATH/perf.dat

    $COLLECT_CMD -o $EXPERIMENT_DIR/blk_${m}_${bs}.er $BIN blk $m $(($m+10)) $(($m+20)) $bs
    er_print -script er_print_script $EXPERIMENT_DIR/blk_${m}_${bs}.er | grep "matmult_[^r-r][^e-e][^f-f]" >> $OUT_PATH/cache.dat
    er_print -script er_print_script $EXPERIMENT_DIR/blk_${m}_${bs}.er | grep "blkmult" >> $OUT_PATH/cache.dat
    truncate -s -1 $OUT_PATH/cache.dat
    echo ",blk,$m,$bs" >> $OUT_PATH/cache.dat
  done
done