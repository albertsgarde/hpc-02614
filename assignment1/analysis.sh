MS="10 50 100 200 500 750 1000 1500"

rm nat_Ofast_interchange.dat
rm lib_Ofast_interchange.dat
rm mnk_Ofast_interchange.dat
rm mkn_Ofast_interchange.dat
rm nmk_Ofast_interchange.dat
rm nkm_Ofast_interchange.dat
rm kmn_Ofast_interchange.dat
rm knm_Ofast_interchange.dat

for m in $MS
do
  ./matmult_c.gcc nat $m $(($m+10)) $(($m+20)) >> nat_Ofast_interchange.dat
  ./matmult_c.gcc lib $m $(($m+10)) $(($m+20)) >> lib_Ofast_interchange.dat
  
  ./matmult_c.gcc mnk $m $(($m+10)) $(($m+20)) >> mnk_Ofast_interchange.dat
  ./matmult_c.gcc mkn $m $(($m+10)) $(($m+20)) >> mkn_Ofast_interchange.dat
  ./matmult_c.gcc nmk $m $(($m+10)) $(($m+20)) >> nmk_Ofast_interchange.dat
  ./matmult_c.gcc nkm $m $(($m+10)) $(($m+20)) >> nkm_Ofast_interchange.dat
  ./matmult_c.gcc kmn $m $(($m+10)) $(($m+20)) >> kmn_Ofast_interchange.dat
  ./matmult_c.gcc knm $m $(($m+10)) $(($m+20)) >> knm_Ofast_interchange.dat
done