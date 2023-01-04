MS="10 50 100 200 500 750 1000 1500"

rm nat_Ofast.dat
rm lib_Ofast.dat
rm mnk_Ofast.dat
rm mkn_Ofast.dat
rm nmk_Ofast.dat
rm nkm_Ofast.dat
rm kmn_Ofast.dat
rm knm_Ofast.dat

for m in $MS
do
  ./matmult_c.gcc nat $m $(($m+10)) $(($m+20)) >> nat_Ofast.dat
  ./matmult_c.gcc lib $m $(($m+10)) $(($m+20)) >> lib_Ofast.dat
  
  ./matmult_c.gcc mnk $m $(($m+10)) $(($m+20)) >> mnk_Ofast.dat
  ./matmult_c.gcc mkn $m $(($m+10)) $(($m+20)) >> mkn_Ofast.dat
  ./matmult_c.gcc nmk $m $(($m+10)) $(($m+20)) >> nmk_Ofast.dat
  ./matmult_c.gcc nkm $m $(($m+10)) $(($m+20)) >> nkm_Ofast.dat
  ./matmult_c.gcc kmn $m $(($m+10)) $(($m+20)) >> kmn_Ofast.dat
  ./matmult_c.gcc knm $m $(($m+10)) $(($m+20)) >> knm_Ofast.dat
done